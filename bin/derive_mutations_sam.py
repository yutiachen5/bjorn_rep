#!/usr/bin/env python

import os
import re
import mappy as mp
import argparse
import polars as pl
import pysam # https://pysam.readthedocs.io/en/latest/api.html


def parse_md_tag(read): 
    # numbers in md is for ref genome
    # case when len(ref) == len(query)
    # pos in the mut file is 1-based

    md = read.get_tag("MD") # 1-based, 22554T: 22555-th nuc of ref is T
    re_expr = r"\d+|\^[A-Z]+|[A-Z]"
    md_tokens = re.findall(re_expr, md)

    ref_pos = read.reference_start # 0-based
    query_pos = read.query_alignment_start # 0-based
    query_seq = read.query_sequence

    token_iter = iter(md_tokens)
    current_token = next(token_iter)

    count  = 0 # record number of mutations
    variants = [] # [(query_pos, ref_seq, query_seq), ...], pos are in 1-bsaed

    for operation, length in read.cigartuples:
        if operation == 0: # alignment
            if current_token.isdigit():
                ref_pos += int(current_token)
                query_pos += int(current_token)
                count += int(current_token)
                current_token = next(token_iter)
            else:
                # deletion
                if current_token.startwith("^"):
                    del_seq = current_token[1:] # skip ^
                    ref_pos += length
                    variants.append((query_pos+1, current_token, del_seq)) # 0-based -> 1-based
                    current_token = next(token_iter)
                # substitution
                else:
                    variants.append((query_pos+1, current_token, query_seq[query_pos])) # 0-based -> 1-based
                    current_token = next(token_iter)

        elif operation == 1: # insertion
            ins_seq = query_seq[query_pos:query_pos+length]
            variants.append(ins_seq)
            query_pos += length
            count += length


        assert count ==  read.get_tag("NM")
    
    print(variants)
    return 0

def parse_sam_file(sam_path):
    sam_file = pysam.AlignmentFile(sam_path, "r")
    unmapped_pos = []
    for read in sam_file:
        # todo: check if the alignment file only contains unmapped reads
        print(read.query_name, read.reference_name,
              read.reference_start, read.reference_end,
              read.query_alignment_start, read.query_alignment_end,read.cigarstring
            )  
        print(read.get_tag("MD"))
        print(read.get_tag("NM"))
        print(read.cigartuples)

        # parse_md_tag(read)
        # parse_variants(read)

    return 0
                

def get_query_mutations(mut, unmapped, seq):
    # for the mapped positions, keep the original mutations in the file
    # for the unmapped positions, compare the nuc in original file with query file
    # if ref.alt_nuc == query.ref_auc: rm the line
    # elif ref.alt_nuc != query.ref_auc: change the nuc_ref from ref.ref_auc to query.ref_nuc

    # create a schema for query genome, pos is 1-based in mutation file
    query_schema = {
        "pos": pl.Int64,
        "query_ref":pl.String
    }

    query_df = pl.DataFrame(
        data={
            "pos": range(1, len(seq)+1),
            "query_ref": list(seq)
        },
        schema = query_schema
    )
    interval_expr = [(pl.col("pos")>=coords[0]) & (pl.col("pos")<coords[1]) for coords in unmapped] # object only, not bool cond
    if interval_expr:
        condition = pl.fold(pl.lit(False), lambda e1, e2: e1 | e2, interval_expr)
        unmatched = mut.filter(condition)
    else:
        unmatched = pl.DataFrame(schema=mut.schema)


    matched = mut.join(unmatched, on=mut.columns, how="anti")
    print(matched.height)
    print(unmatched.height)
    print(mut.height)

    unmatched_query_ref_nuc = unmatched.join(query_df, on='pos')
    unmatched_query_ref_nuc = unmatched_query_ref_nuc.with_columns(
        pl.when(pl.col("alt") == pl.col("query_ref"))
            .then(pl.lit("-999"))
            .otherwise(pl.col("query_ref"))
            .alias("alt")
    )
    unmatched = unmatched_query_ref_nuc.filter(
        pl.col("alt") != "-999"
    ).drop("query_ref")
    print(unmatched.filter((pl.col('pos') == 22578) & (pl.col('sra') == 'hCoV-19/USA/STM-E4PG4VV42/2021')))

    query_mut = pl.concat([matched, unmatched], how='vertical').sort("pos")

    print(query_mut.height)

    return query_mut



def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("--mut", help="Input mutation file path in TSV format.", required=True)
    # parser.add_argument("--ref", help="Reference genome which derived the mutation file.", required=True)
    # parser.add_argument("--query", help="Query genome.", required=True)
    parser.add_argument('-s', help="SAM file with tags for the alignment between reference and query genome.", required=True)
    parser.add_argument("-o", help="Output diractory for new mutation file.")
    args = parser.parse_args()

    parse_sam_file(args.s)

    # unmapped, record, query_seq = alignment(args.ref, args.query)

    # mut = pl.read_csv(args.mut, separator="\t", has_header=True)

    # query_mut = get_query_mutations(mut, unmapped, query_seq)

    # query_mut.write_csv(os.path.join(args.output_dir, record + "_mutations.tsv"), separator="\t")

    
if __name__ == "__main__":
    main()
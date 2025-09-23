#!/usr/bin/env python

import argparse
import polars as pl
from Bio import SeqIO

def get_seq_from_reference(reference: str, pos: int, length: int) -> str:
    """
    pos is 1-based
    """
    if pos is None:
        return None, None
    pos -= 1
    if reference is None:
        return None, None
    if pos + length > len(reference):
        return None, None
    return reference[pos], reference[pos:pos + length]

def parse_mutations(df: pl.DataFrame, mutation_col: str, region: str = None, sra_regex: str=  r"(.*)", reference_fasta: str= None) -> pl.DataFrame:
    """
    Parse mutation strings into structured columns using Polars native expressions.

    The following patterns for mutations are expected:

    nuc:A27259C
    aa:YP_009724394.1:D61L(nuc:G27382C;nuc:A27383T;nuc:T27384C)(ref:GAT-alt:CTC)
    aa:YP_009724397.2:P13L(nuc:C28311T)(ref:CCC-alt:CTC)
    del:28362:9
    ins:2028:3

    Args:
        df: Input dataframe with mutation column
        mutation_col: Name of column containing mutation strings
        region: Value to fill in region column (default None)
    
    Returns:
        DataFrame with parsed mutation columns
    """
    # Explode mutations
    df = df.with_columns(
        pl.col("mutations").str.split("|")
    ).explode("mutations") # expand each element of the list into its own row

    # Define regex patterns
    PATTERNS = {
        'mut_type': r"^([^:]+):",
        'base_change': r"[A-Z](\d+)[A-Z]",  # captures position from nucleotide change
        'nuc_change': r"([A-Z])(\d+)([A-Z])",  # captures ref, pos, alt from nucleotide
        'position_after': r":(\d+):",  # captures position after colon
        'gff_feature': r"aa:([^:]+):",
        'aa_change': r"aa:[^:]+:([A-Z])(\d+)([A-Z])",  # captures aa ref, pos, alt
        'codon_pair': r"ref:([A-Z]{3})-alt:([A-Z]{3})",  # captures both codons
        'nuc_extract_all': r"nuc:([A-Z]\d+[A-Z])" # Captures individual nucleotide changes
    }
    
    # Define empty dataframe schema once
    EMPTY_SCHEMA = {
        "sra": pl.String, "region": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String,
        "GFF_FEATURE": pl.String, "ref_codon": pl.String, "alt_codon": pl.String,
        "ref_aa": pl.String, "alt_aa": pl.String, "pos_aa": pl.Int64
    }
    
    # Step 1: Extract mutation type (first part before colon)
    result = df.with_columns([
        pl.col(mutation_col).str.extract(PATTERNS['mut_type']).alias("mut_type"),
        pl.lit(region).alias("region"), # add a new column
        pl.col("query").str.extract(sra_regex).alias("sra")
    ])
    
    # Step 2: Filter and process each mutation type separately to avoid creating intermediate columns
    
    # Extract ALL nucleotide changes from aa entries for later processing
    result = result.with_columns([
        pl.when(pl.col("mut_type") == "aa").then(
            pl.col(mutation_col).str.extract_all(PATTERNS['nuc_extract_all'])
        ).alias("nuc_matches")
    ]).with_columns([
        pl.col("nuc_matches").list.len().fill_null(0).alias("nuc_count")
    ])
    
    # Step 3: Process ALL aa entries (single and multiple nucleotides) the same way
    aa_entries = result.filter(
        (pl.col("mut_type") == "aa") & 
        (pl.col("nuc_count") > 0)
    )
    
    if len(aa_entries) > 0:
        # Explode all aa_entries
        aa_final = aa_entries.with_columns([
            pl.col("nuc_matches").list.eval(
                pl.element().str.extract_groups(PATTERNS['nuc_change'])
            ).alias("nuc_components"),
            pl.col(mutation_col).str.extract(PATTERNS['gff_feature']).alias("GFF_FEATURE"),
            pl.col(mutation_col).str.extract_groups(PATTERNS['codon_pair']).alias("codon_parts"),
            pl.col(mutation_col).str.extract_groups(PATTERNS['aa_change']).alias("aa_parts")
        ]).with_columns([
            pl.col("nuc_components").list.eval(pl.element().struct.field("1")).alias("aa_nuc_refs"),
            pl.col("nuc_components").list.eval(pl.element().struct.field("2").cast(pl.Int64)).alias("aa_nuc_positions"),
            pl.col("nuc_components").list.eval(pl.element().struct.field("3")).alias("aa_nuc_alts")
        ]).explode([
            "aa_nuc_refs", "aa_nuc_positions", "aa_nuc_alts"
        ]).with_columns([
            pl.col("region").alias("region"),
            pl.col("aa_nuc_positions").alias("pos"),
            pl.col("aa_nuc_refs").alias("ref"),
            pl.col("aa_nuc_alts").alias("alt"),
            pl.col("GFF_FEATURE"),
            pl.col("codon_parts").struct.field("1").alias("ref_codon"),
            pl.col("codon_parts").struct.field("2").alias("alt_codon"),
            pl.col("aa_parts").struct.field("1").alias("ref_aa"),
            pl.col("aa_parts").struct.field("2").cast(pl.Int64).alias("pos_aa"),
            pl.col("aa_parts").struct.field("3").alias("alt_aa")
        ]).select(list(EMPTY_SCHEMA.keys()))
    else:
        aa_final = pl.DataFrame(schema=EMPTY_SCHEMA)
    
    # Step 4: Handle simple cases (nuc, del, ins)
    simple_cases = result.filter(pl.col("mut_type").is_in(["nuc", "del", "ins"]))

    
    if len(simple_cases) > 0:
        EMPTY_AA_FIELDS = ["GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa"]

        simple_final = simple_cases.with_columns([            
            # Extract all nuc components at once, then handle position for del/ins
            pl
            .when(pl.col("mut_type") == "nuc")
            .then(pl.col(mutation_col).str.extract_groups(f"nuc:{PATTERNS['nuc_change']}")
                    .struct.rename_fields(["ref", "pos", "alt"])
                    .struct.with_fields([
                            pl.field("pos").cast(pl.Int64),
                            pl.field("ref"),
                            pl.field("alt")
                        ])
                  )
            .when(pl.col("mut_type") == "del")
            .then(
                pl.col(mutation_col).str.extract_groups(r'^[^:]+:(\d+):(\d+)$')
                .alias("tmp_extract")
                .struct.rename_fields(["pos", "length"])
                .struct.with_fields([
                    pl.field("pos").cast(pl.Int64),
                    pl.field("length").cast(pl.Int64)
                ])
                .map_elements(
                    lambda x: (lambda del_ref_result ,seq_result: {
                        "pos": x["pos"] - 1 if x["pos"] is not None else None, #TODO: This is hacky, but works for now
                        "ref": del_ref_result[0],
                        "alt": "-" + seq_result[1] if seq_result[1] else None
                    })(get_seq_from_reference(reference_fasta, x["pos"] - 1 if x["pos"] is not None else None, 1), get_seq_from_reference(reference_fasta, x["pos"], x["length"])),
                    return_dtype=pl.Struct([
                        pl.Field("pos", pl.Int64),
                        pl.Field("ref", pl.String),
                        pl.Field("alt", pl.String)
                    ]))
                )
            .when(pl.col("mut_type") == "ins")# TODO: Figure out how to get inserted bases for insertions
            .then(
                pl.col(mutation_col).str.extract_groups(r'^[^:]+:(\d+):(\d+)$')
                .alias("tmp_extract")
                .struct.rename_fields(["pos", "length"])
                .struct.with_fields([
                    pl.field("pos").cast(pl.Int64),
                    pl.field("length").cast(pl.Int64)
                ])
                .map_elements(
                    lambda x: (lambda seq_result: {
                        "pos": x["pos"],
                        "ref": seq_result[0], 
                        "alt": f"+{x['length']}"
                    })(get_seq_from_reference(reference_fasta, x["pos"], x["length"])),
                    return_dtype=pl.Struct([
                        pl.Field("pos", pl.Int64),
                        pl.Field("ref", pl.String),
                        pl.Field("alt", pl.String)
                    ]))
                )
            .alias("components")
        ]).with_columns([
            # Extract from struct
            pl.col("components").struct.field("pos").cast(pl.Int64).alias("pos"),
            pl.col("components").struct.field("ref").alias("ref"),
            pl.col("components").struct.field("alt").alias("alt"),
            
            # Empty strings for other columns
            *[pl.lit(None, dtype=pl.String).alias(col) for col in EMPTY_AA_FIELDS],
            pl.lit(None, dtype=pl.Int64).alias("pos_aa")
        ]).select(list(EMPTY_SCHEMA.keys()))
        
    else:
        simple_final = pl.DataFrame(schema=EMPTY_SCHEMA)
    
    # Step 5: Combine all results
    final_result = pl.concat([simple_final, aa_final], how="vertical")
    
    return final_result

def read_reference(fasta_path: str) -> list[str, str]:
    """
    Read a reference FASTA file and return the sequence as a string.
    
    Args:
        fasta_path: Path to the reference FASTA file.
    
    Returns:
        The concatenated sequence from the FASTA file.
    """
    try:
        fasta_iter = SeqIO.parse(fasta_path, "fasta")
        try:
            record = next(fasta_iter)
        except StopIteration:
            e = f"Error: No records found in FASTA file: {fasta_path}"
            return None, e
        return str(record.seq), None
    except Exception as e:
        return None, f"Error reading FASTA file {fasta_path}"

def main():
    parser = argparse.ArgumentParser(description="Convert gofasta variant CSVs to a single TSV.")
    parser.add_argument("--csv_path", help="Input CSV file path", required=True)
    parser.add_argument("--output", default="mutations.tsv", help="Output TSV file name.")
    parser.add_argument("--region", default=None, help="Region name", required=True)
    parser.add_argument("--reference", default=r"(.*)", help="Reference fasta. Needed to populate insertions and deletions.", required=True)
    parser.add_argument("--sra-regex", default=r"(.*)", help="Regex to get SRA")
    args = parser.parse_args()

    df = pl.read_csv(args.csv_path)

    if df.is_empty():
        exit(f"Error: {args.csv_path} is empty.")
                
    if len(df.columns) != 2:
        exit(f"Warning: Expected 2 columns in {args.csv_path}, found {len(df.columns)}")

    reference_fasta, e = read_reference(args.reference)
    if e:
        exit(e)

    result = parse_mutations(df, "mutations", args.region, args.sra_regex, reference_fasta)
    result.sort("pos").write_csv(
        args.output,
        separator="\t",
        include_header=True
    )

    
if __name__ == "__main__":
    main()
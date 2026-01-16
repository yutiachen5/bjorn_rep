from Bio import SeqIO
import os

# COORDS asre base 1
SPIKE_COORDS = [21563, 25384] # coords are from NC_04512.2 gff file for spike protein
RBD_COORDS = [22553, 23155] # 22553 = 21563 + (331 - 1) * 3

def apply_dms_ref(ref, gene_seq, out_name, output_prefix, region = "spike"):
    region_coords = SPIKE_COORDS if region == "spike" else RBD_COORDS
    hu1_with_dms_ref = ref.seq[:region_coords[0]-2] + gene_seq + ref.seq[region_coords[1]:]
    hu1_with_dms_ref.id = out_name
    hu1_with_dms_ref.name = out_name
    hu1_with_dms_ref.description = out_name
    SeqIO.write(hu1_with_dms_ref, os.path.join(output_prefix, out_name) , format="fasta")

def get_gene_from_dms_ref(dms_ref, feature_type="gene"):
    gene_feature = next(i for i in dms_ref.features if i.type == feature_type)
    return dms_ref[gene_feature.location.start: gene_feature.location.end]

ref = SeqIO.read("./data/hu1/NC_045512.2.fasta", format="fasta")

#path_to_dms_ref = "./data/xbb.1.5/rbd/PacBio_amplicon.gb"
#path_to_dms_ref = "./data/xbb.1.5/spike/PacBio_amplicon.gb"
path_to_dms_ref = "./data/ba.2/spike/PacBio_amplicon.gb"
dms_ref = SeqIO.read(path_to_dms_ref, format="genbank")

gene_seq = get_gene_from_dms_ref(dms_ref)
apply_dms_ref(ref, gene_seq, "NC_045512.2_BA.2_spike.fa", "./output/ba.2/spike")
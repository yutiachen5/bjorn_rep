from Bio import AlignIO, SeqIO
import os

# 1-based coords
SPIKE_COORDS = [21563, 25384]
RBD_COORDS = [22553, 23155] # 22553 = 21563 + 330 *3

def apply_dms_ref(ref, gene_seq, out_name, output_prefix, region = "spike"):
    region_coords = SPIKE_COORDS if region == "spike" else RBD_COORDS
    hu1_with_dms_ref = ref.seq[:region_coords[0]-1] + gene_seq + ref.seq[region_coords[1]:]
    hu1_with_dms_ref.id = out_name
    hu1_with_dms_ref.name = out_name
    hu1_with_dms_ref.description = out_name
    SeqIO.write(hu1_with_dms_ref, os.path.join(output_prefix, out_name) , format="fasta")


ref = SeqIO.read("./data/hu1/NC_045512.2.fasta", format="fasta")

# aligned by mafft
# algn = AlignIO.read("./data/BA.1-BA.2-BA.4-BA.5-BA.2.75.fa", format="fasta") 
algn = AlignIO.read("./data/JN.1-KP.2-KP.3-XBB.1.5.fa", format="fasta") 



# coords shift because of ins/del in alignemnt
# rbd_algn = algn[:, 22561: 23164] 
rbd_algn = algn[:, 22564: 23167] 


for rbd in rbd_algn[1:]:
    out_name = "NC_045512.2_escape_{}_rbd.fa".format(rbd.id.split("|")[0])
    apply_dms_ref(ref, rbd, out_name, "./output/", region = "rbd")
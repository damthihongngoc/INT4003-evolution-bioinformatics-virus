from Bio import SeqIO

REF_ID = "MN975262.1"
FA = "data/genomic_aln.fa"

#  Định nghĩa gen và tọa độ gen, lấy tọa độ gene trên genome tham chiếu
genes = {
    "ORF1ab": (266, 21555),
    "S": (21563, 25384),
    "ORF3a": (25393, 26220),
    "E": (26245, 26472),
    "M": (26523, 27191),
    "ORF6": (27202, 27387),
    "ORF7a": (27394, 27759),
    "ORF8": (27894, 28259),
    "N": (28274, 29533),
    "ORF10": (29558, 29674),
}

# Đọc alignment và lấy sequence tham chiếu 
records = SeqIO.to_dict(SeqIO.parse(FA, "fasta"))
seq = records[REF_ID].seq

genome_pos = 0
mapping = {}
#Duyệt từng cột alignment, nếu không phải gap, tăng vị trí genome
for aln_pos, base in enumerate(seq, start=1):
    if base != "-":
        genome_pos += 1
        mapping[genome_pos] = aln_pos

with open("data/genomic_aln.nex", "w") as out:
    out.write("#nexus\n")
    out.write("begin sets;\n")
    for gene, (start, end) in genes.items():
        out.write(
            f"  charset {gene} = {mapping[start]}-{mapping[end]};\n"
        )
    out.write("end;\n")

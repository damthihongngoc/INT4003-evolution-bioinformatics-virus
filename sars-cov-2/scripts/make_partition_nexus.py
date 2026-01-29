from Bio import SeqIO

REF_ID = "MN975262.1"
FA = "data/genomic_aln.fa"

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

records = SeqIO.to_dict(SeqIO.parse(FA, "fasta"))
seq = records[REF_ID].seq

genome_pos = 0
mapping = {}

for aln_pos, base in enumerate(seq, start=1):
    if base != "-":
        genome_pos += 1
        mapping[genome_pos] = aln_pos

with open("results/sarscov2_partition.nex", "w") as out:
    out.write("#nexus\n")
    out.write("begin sets;\n")
    for gene, (start, end) in genes.items():
        out.write(
            f"  charset {gene} = {mapping[start]}-{mapping[end]};\n"
        )
    out.write("end;\n")

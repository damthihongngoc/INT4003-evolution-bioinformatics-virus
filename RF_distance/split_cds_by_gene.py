from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os

core_genes = {"ORF1ab", "S", "E", "M", "N"}
genes = defaultdict(list)

for record in SeqIO.parse("sars2_26.gb", "genbank"):
    acc = record.id

    for feat in record.features:
        if feat.type != "CDS":
            continue

        gene = feat.qualifiers.get("gene", [""])[0]
        if gene not in core_genes:
            continue

        seq = feat.extract(record.seq)

        seqrec = SeqRecord(
            seq,
            id=f"{acc}|{gene}",
            description=""
        )

        genes[gene].append(seqrec)

os.makedirs("genes", exist_ok=True)

for g, seqs in genes.items():
    SeqIO.write(seqs, f"genes/{g}.fasta", "fasta")

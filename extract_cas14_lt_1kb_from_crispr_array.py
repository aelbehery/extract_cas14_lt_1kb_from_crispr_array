#!/usr/bin/env python3

import sys
import csv
from Bio import SeqIO
import gzip

contigs = []
correct_cas14 = []
cas14_seqs = []
correct_cas14_coords = []
contig_desc = []
assembly_accession = "_".join(sys.argv[1].split("/")[1].split("_")[0:2])


with open(sys.argv[1]) as gff:
    reader = list(csv.reader(gff, delimiter='\t', quoting=csv.QUOTE_NONE))
    for row in reader:
        if not row[0].startswith("#"):
            contigs.append(row[0])
contigs = set(contigs)

for contig in contigs:
    cas14 = []
    dr = []
    for row in reader:
        if row[0] == contig:
            if "cas14" in row[8]:
                cas14.append([row[3],row[4],row[6]])
            elif row[2] == "direct_repeat":
                dr.append([row[3],row[4],row[8].split(";")[1].split(":")[1]])
    for protein in cas14:
        if protein[2] == "+":
            diff = [int(x[0]) - int(protein[1]) for x in dr]
        else:
            diff = [int(protein[0]) - int(x[1]) for x in dr]
        positions = [index for index, element in enumerate(diff) if 0 < element <= 1000]
        if len(positions) == 1:
            position = positions[0]
            correct_cas14.append([contig] + protein + dr[position] + [str(diff[position])])


if correct_cas14:
    for i in range(len(correct_cas14)):
        correct_cas14_coords.append("-".join(correct_cas14[i][1:3]))
    prots = list(SeqIO.parse(sys.argv[1].rsplit(".",1)[0]+".faa", "fasta"))
    for correct_cas14_coord in correct_cas14_coords:
        for prot in prots:
            if correct_cas14_coord in prot.description:
                cas14_seqs.append(prot.seq)
    with gzip.open(sys.argv[1].split("/")[1]+"_genomic.fna.gz", "rt") as contigs_fasta:
        contig_records = list(SeqIO.parse(contigs_fasta, "fasta"))
    for contig in contigs:
        for contig_record in contig_records:
            if contig == contig_record.id:
                contig_desc.append(contig_record.description.split(" ", 1)[1])
    for correct_cas14_prot in correct_cas14:
        j = 0
        print(assembly_accession, contig_desc[j], "\t".join(correct_cas14_prot), cas14_seqs[j], len(cas14_seqs[j]), sep = "\t")
        j +=1

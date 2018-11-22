#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description='Genotype VGSC gene.')

parser.add_argument("ref_path",
                   help = 'Reference genome fasta path')

parser.add_argument("mdom_path",
                   help = 'M. domestica aa fasta path')

parser.add_argument("bed_path",
                   help = "Reference bed path")

args = parser.parse_args()

muscle_path = "muscle"
ref_path = args.ref_path
bed_path = args.bed_path
Mdom_path = args.mdom_path

dna = {seq.id:seq.seq for seq in SeqIO.parse(ref_path, "fasta")}

with open(bed_path) as f:
    bed = [l.rstrip().split("\t")  for l in f.readlines()]

vgsc_ck = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not l[3] in ["Exon19d", "Exon26l"]]
vgsc_dl = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not l[3] in ["Exon19c", "Exon26k"]]

if bed[0][3] == "Exon1":
    AA_ck = str(sum(vgsc_ck, Seq("")).translate())
    AA_dl = str(sum(vgsc_dl, Seq("")).translate())
else:
    AA_ck = str(sum(vgsc_ck, Seq("")).reverse_complement().translate())
    AA_dl = str(sum(vgsc_dl, Seq("")).reverse_complement().translate())

Mdom_AA = [seq.seq for seq in SeqIO.parse(Mdom_path, "fasta")]
Mdom_AA = str(Mdom_AA[0])

proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n>Mos\n" + AA_ck], stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path], stdin = proc1.stdout, stdout = subprocess.PIPE)

aligned = [line.decode().rstrip() for line in proc2.stdout.readlines()]

for l in aligned:
    print(l)

sys.exit(0)

Mdom_AA_aligned = []
Mos_AA_aligned = []
is_Mdom = True
for line in aligned[1:]:
    if line.startswith(">"):
        is_Mdom = False
        continue
    if is_Mdom:
        Mdom_AA_aligned += list(line)
    else:
        Mos_AA_aligned += list(line)

if not len(Mdom_AA_aligned) == len(Mos_AA_aligned):
    raise ValueError("error!")

idx = [0, 0]

for mos, mdom in zip(Mos_AA_aligned, Mdom_AA_aligned):
    if mos != "-":
        idx[0] += 1
        a = str(idx[0])
    else:
        a = "-"
    if mdom != "-":
        idx[1] += 1
        b = str(idx[1])
    else:
        b = "-"

    print(a + "\t.\t "+ b)

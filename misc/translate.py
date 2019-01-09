#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import subprocess
import argparse
import sys
import json

parser = argparse.ArgumentParser(description='Genotype VGSC gene.')

parser.add_argument("ref_path",
                   help = 'Reference genome fasta path')

parser.add_argument("bed_path",
                   help = "Reference bed path")

parser.add_argument("out_fasta_path",
                   help = "fasta to output")

parser.add_argument("-m", "--mdom_path", dest = "mdom_path",
                   help = 'M. domestica aa fasta path',
                   default = sys.path[0] + "/AAB47604.fa")

args = parser.parse_args()

muscle_path = "muscle"
ref_path = args.ref_path
bed_path = args.bed_path
Mdom_path = args.mdom_path

dna = {seq.id:seq.seq for seq in SeqIO.parse(ref_path, "fasta")}

with open(bed_path) as f:
    bed = [l.rstrip().split("\t")  for l in f.readlines()]

vgsc_ck = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not (l[3].endswith('d') or  l[3].endswith('l'))]
vgsc_dl = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not (l[3].endswith('c') or  l[3].endswith('k'))]

if bed[0][3] == "Exon1":
    AA_ck = str(sum(vgsc_ck, Seq("")).translate())
    AA_dl = str(sum(vgsc_dl, Seq("")).translate())
else:
    AA_ck = str(sum(vgsc_ck, Seq("")).reverse_complement().translate())
    AA_dl = str(sum(vgsc_dl, Seq("")).reverse_complement().translate())

Mdom_AA = [seq.seq for seq in SeqIO.parse(Mdom_path, "fasta")]
Mdom_AA = str(Mdom_AA[0])

aligned = [ "", ""]
proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n>Mos\n" + AA_ck], stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path], stdin = proc1.stdout, stdout = subprocess.PIPE)

aligned[0] = [line.decode().rstrip() for line in proc2.stdout.readlines()]

proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n>Mos\n" + AA_dl], stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path], stdin = proc1.stdout, stdout = subprocess.PIPE)

aligned[1] = [line.decode().rstrip() for line in proc2.stdout.readlines()]

with open(args.out_fasta_path, 'w') as f:
    for l in aligned[0]:
        print(l, file = f)

for l in aligned[1]:
    print(l, file = sys.stderr)

for i in (0 , 1):

    Mdom_AA_aligned = []
    Mos_AA_aligned = []
    is_Mdom = True
    for line in aligned[i][1:]:
        if line.startswith(">"):
            is_Mdom = False
            continue
        if is_Mdom:
            Mdom_AA_aligned += list(line.rstrip())
        else:
            Mos_AA_aligned += list(line.rstrip())

    if not len(Mdom_AA_aligned) == len(Mos_AA_aligned):
        raise ValueError("error!")

    idx = [0, 0]

    with open(sys.path[0] + "/../scripts/kdr_list.json") as f:
        kdr_dict = json.load(f)

    for mos, mdom in zip(Mos_AA_aligned, Mdom_AA_aligned):
        kdr_symbol = ""
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

        if mos != "-" and mdom != "-":
            if mos != mdom:
                if (b + mos) in kdr_dict:
                    if kdr_dict[(b + mos)] == 1:
                        kdr_symbol = "kdr!"
                    else:
                        kdr_symbol = "kdr?"
                print(b + "\t" + mdom + "->" + mos + "  " + kdr_symbol, file = sys.stderr)

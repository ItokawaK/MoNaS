#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import subprocess
import argparse
import sys
import json
import os
from configuration import GenomeRef

def show_alignment(algn, file):

    seq_list =[[], []]
    names = ["", ""]
    idx = -1
    for l in algn:
        if l.startswith(">"):
            idx += 1
            names[idx] = l.lstrip(">").rstrip()
        else:
            seq_list[idx].append(l.rstrip())

    for i in range(len(seq_list[0])):
        seq1 = seq_list[0][i]
        seq2 = seq_list[1][i]
        matches = ""
        for ii in range(len(seq1)):
            match = " "
            if seq1[ii] is not "-" and seq2 is not "-":
                if seq1[ii] == seq2[ii]:
                    match = "|"
            matches += match

        print(names[0] + "\t" + seq1, file = file)
        print("    " + "\t" + matches, file = file)
        print(names[1] + "\t" + seq2, file = file)
        print("", file = file)




script_dir = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description='Genotype VGSC gene.')

parser.add_argument("-r", "--ref_dir", dest = "ref_dir",
                   help = 'Reference genome roor_dir',
                   default = script_dir + "/../references")

parser.add_argument("-s", "--sepcies", dest = "species",
                   help = "species")

parser.add_argument("-o", "--out_fasta", dest = "out_fasta_path",
                   help = "fasta to output [species/ref.modom.fa]")

parser.add_argument("-m", "--mdom_path", dest = "mdom_path",
                   help = 'M. domestica aa fasta path',
                   default = script_dir + "/../misc/AAB47604.fa")

args = parser.parse_args()

gref = GenomeRef(args.ref_dir, args.species)

muscle_path = "muscle"
ref_path = gref.ref_fa
bed_path = gref.bed
Mdom_path = args.mdom_path

if not args.out_fasta_path:
    out_fasta = os.path.join(args.ref_dir, args.species, "ref.mdom.fa")
else:
    out_fasta = args.out_fasta_path

dna = {seq.id:seq.seq for seq in SeqIO.parse(ref_path, "fasta")}

with open(bed_path) as f:
    bed = [l.rstrip().split("\t")  for l in f.readlines()]

vgsc_ck = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not (l[3].endswith('d') or  l[3].endswith('l'))]
vgsc_dl = [dna[l[0]][int(l[1]):int(l[2])] for l in bed if not (l[3].endswith('c') or  l[3].endswith('k'))]

if bed[0][3] == "Exon1":
    AA_ck = str(sum(vgsc_ck, Seq("")).translate()).rstrip("*")
    AA_dl = str(sum(vgsc_dl, Seq("")).translate()).rstrip("*")
else:
    AA_ck = str(sum(vgsc_ck, Seq("")).reverse_complement().translate()).rstrip("*")
    AA_dl = str(sum(vgsc_dl, Seq("")).reverse_complement().translate()).rstrip("*")

cnt = 0
for aa in AA_ck:
    cnt +=1
    if aa is "*":
        print("Warning !: Imature stop codon (*) was found at the " + cnt +
              "-th position in in VGSC_ck.", file = sys.stderr)
cnt = 0
for aa in AA_dl:
    cnt +=1
    if aa is "*":
        print("Warning !: Imature stop codon (*) was found at the " + cnt +
              "-th position in in VGSC_dl.", file = sys.stderr)


Mdom_AA = [seq.seq for seq in SeqIO.parse(Mdom_path, "fasta")]
Mdom_AA = str(Mdom_AA[0])

#show_alignment(aligned[0], sys.stderr)

proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + \
                                  "\n>Mos_ck\n" + AA_ck + \
                                  "\n>Mos_dl\n" + AA_dl],
                         stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path, "-clw"], stdin = proc1.stdout, stdout = subprocess.PIPE)

for l in proc2.stdout.readlines():
    print(l.decode().rstrip(), file=sys.stderr)

print("Do you accept this alignment? 0: no, >0: yes", file = sys.stderr)
to_proceed = input(">>>>  ")
if to_proceed == "0":
    sys.exit(1)


aligned = [ "", ""]
proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n>Mos\n" + AA_ck], stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path], stdin = proc1.stdout, stdout = subprocess.PIPE)

aligned[0] = [line.decode().rstrip() for line in proc2.stdout.readlines()]

proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n>Mos\n" + AA_dl], stdout = subprocess.PIPE)
proc2 = subprocess.Popen([muscle_path], stdin = proc1.stdout, stdout = subprocess.PIPE)

aligned[1] = [line.decode().rstrip() for line in proc2.stdout.readlines()]

with open(out_fasta, 'w') as f:
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

    kdr_list = os.path.dirname(os.path.abspath(__file__)) + "/../scripts/kdr_list.json"

    with open(kdr_list) as f:
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

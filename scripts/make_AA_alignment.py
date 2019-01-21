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
import extract_exons

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description='Genotype VGSC gene.')

    parser.add_argument("-r", "--ref_dir", dest = "ref_dir",
                       help = 'Reference genome roor_dir [MoNaS/references] ',
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

    for file in [ref_path, bed_path, Mdom_path]:
        if not os.path.isfile(file):
            print(file + " was not found!", file = sys.stderr)
            sys.exit(1)

    if not args.out_fasta_path:
        out_fasta = os.path.join(args.ref_dir, args.species, "ref.mdom.fa")
    else:
        out_fasta = args.out_fasta_path

    dna = extract_exons.extract(ref_path, bed_path)

    vgsc_mrna = dict()
    vgsc_prot = dict()

    vgsc_mrna['ck'] = [exon for exon in dna if not exon.id[-1] in  ["d", "l"]]
    vgsc_mrna['dl'] = [exon for exon in dna if not exon.id[-1] in  ["c", "k"]]

    for variant in ['ck', 'dl']:
        vgsc_prot[variant] = str(sum(vgsc_mrna[variant], Seq("")).translate().seq).rstrip("*")

        cnt = 0
        for aa in vgsc_prot[variant]:
            cnt += 1
            if aa is "*":
                print("Warning !: Imature stop codon (*) was found at the "
                      + str(cnt)
                      + "-th position in in VGSC_"
                      + variant + ".",
                      file = sys.stderr)

    Mdom_AA = [seq.seq for seq in SeqIO.parse(Mdom_path, "fasta")]
    Mdom_AA = str(Mdom_AA[0])

    # proc1 = subprocess.Popen(["echo", ">Mdom\n" + Mdom_AA + "\n"\
    #                                   ">Mos_ck\n" + vgsc_prot['ck'] + "\n"\
    #                                   ">Mos_dl\n" + vgsc_prot['dl']
    #                           ],
    #                           stdout = subprocess.PIPE)

    proc = subprocess.Popen([muscle_path, "-clw"],
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE)

    in_fasta = ">Mdom\n" + Mdom_AA + "\n"\
               ">Mos_ck\n" + vgsc_prot['ck'] + "\n"\
               ">Mos_dl\n" + vgsc_prot['dl']

    proc.stdin.write(in_fasta.encode())
    proc.stdin.close()

    for l in proc.stdout.readlines():
        print(l.decode().rstrip(), file=sys.stderr)

    print("Do you accept this alignment? 0: no, >0: yes", file = sys.stderr)
    to_proceed = input(">>>>  ")
    if to_proceed == "0":
        sys.exit(1)

    proc = subprocess.Popen([muscle_path,
                             "-out",  out_fasta],
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)

    proc.stdin.write(in_fasta.encode())
    proc.stdin.close()

    alignment = [seq for seq in SeqIO.parse(out_fasta, 'fasta')]

    kdr_list = os.path.dirname(os.path.abspath(__file__)) + "/../scripts/kdr_list.json"
    with open(kdr_list) as f:
        kdr_dict = json.load(f)

    indecies = [0, 0, 0]
    for pos in range(len(alignment[0].seq)):
        for seq in [0, 1, 2]:
            if alignment[seq][pos] is not "-":
                indecies[seq] += 1
        AA_Change_ck = ""
        AA_Change_dl = ""
        is_kdr_ck = ""
        is_kdr_dl = ""
        if alignment[0][pos] != alignment[1][pos] and\
                    not (alignment[0][pos] == "-" or alignment[1][pos] == "-"):
            AA_Change_ck = "{}{}".format(indecies[0], alignment[1][pos])
            if AA_Change_ck in kdr_dict:
                if kdr_dict[AA_Change_ck] == 1:
                    is_kdr_ck = "kdr!"
                else:
                    is_kdr_ck = "kdr?"
        if alignment[0][pos] != alignment[2][pos] and\
                    not (alignment[0][pos] == "-" or alignment[2][pos] == "-"):
            AA_Change_dl = "{}{}".format(indecies[0], alignment[2][pos])
            if AA_Change_dl in kdr_dict:
                if kdr_dict[AA_Change_dl] == 1:
                    is_kdr_dl = "kdr!"
                else:
                    is_kdr_dl = "kdr?"

        if AA_Change_ck and AA_Change_ck == AA_Change_dl:
            print("{} in ck and dl {}".format(AA_Change_ck, is_kdr_ck), file=sys.stderr)
        elif AA_Change_ck:
            print("{} in ck {}".format(AA_Change_ck, is_kdr_ck), file=sys.stderr)
        elif AA_Change_dl:
            print("{} in dl {}".format(AA_Change_dl, is_kdr_dl), file=sys.stderr)

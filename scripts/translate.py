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


    aligned = [ "", ""]

    for i in (0, 1):
        in_fasta =  ">Mdom\n" + Mdom_AA + "\n"\
                    ">Mos\n" + vgsc_prot[['ck', 'dl'][i]]

        proc = subprocess.Popen([muscle_path],
                                 stdin = subprocess.PIPE,
                                 stdout = subprocess.PIPE)

        proc.stdin.write(in_fasta.encode())
        proc.stdin.close()

        aligned[i] = [line.decode().rstrip() for line in proc.stdout.readlines()]

    with open(out_fasta, 'w') as f:
        for l in aligned[0]:
            print(l, file = f)

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

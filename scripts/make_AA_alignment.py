#!/usr/bin/env python
'''
Copyright (c) 2019, Kentaro Itokawa <itokawa@nih.go.jp>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import subprocess
import argparse
import sys
import json
import os
import extract_exons

def main(ref_path, bed_path, Mdom_path, out_fasta_path = None, translate = None):

    muscle_path = "muscle"

    #Extract CDS sequences, join into mRNA, and translate

    dna = extract_exons.extract(ref_path, bed_path)

    vgsc_mrna = dict()
    vgsc_prot = dict()

    vgsc_mrna['ck'] = [exon for exon in dna if not exon.id[-1] in  ["d", "l"]]
    vgsc_mrna['dl'] = [exon for exon in dna if not exon.id[-1] in  ["c", "k"]]

    for variant in ['ck', 'dl']:
        vgsc_prot[variant] = sum(vgsc_mrna[variant], Seq("")).translate()
        vgsc_prot[variant].seq = vgsc_prot[variant].seq.rstrip("*")
        vgsc_prot[variant].id = "VGSC_" + variant
        vgsc_prot[variant].description = ""

        cnt = 0
        for aa in vgsc_prot[variant]:
            cnt += 1
            if aa != "*":
                print("Warning !: Imature stop codon (*) was found at the "
                      + str(cnt)
                      + "-th position in in VGSC_"
                      + variant + ".",
                      file = sys.stderr)

    vgsc_prot = [vgsc_prot[v] for v in ['ck', 'dl']] # dict to list



    if translate:
        with open(translate, 'w') as f:
            SeqIO.write(vgsc_prot, f, 'fasta')

    #Conduct multiple aligment

    if not out_fasta_path:
        sys.exit(0)

    info = open(out_fasta_path + ".info", 'w')

    Mdom_AA = [seq for seq in SeqIO.parse(Mdom_path, "fasta")]
    vgsc_prot = Mdom_AA + vgsc_prot

    proc = subprocess.Popen([muscle_path, "-clw"],
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE,
                             universal_newlines=True)

    SeqIO.write(vgsc_prot, proc.stdin, "fasta")
    proc.stdin.close()

    print("1. Multiple aligment in ClustW format", file=info)
    for l in proc.stdout.readlines():
        print("    " + l.rstrip(), file=info)
    print("\n2. AAs different from M. domestica VGSC", file=info)

    proc = subprocess.Popen([muscle_path],
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE,
                             universal_newlines=True)

    SeqIO.write(vgsc_prot, proc.stdin, "fasta")
    proc.stdin.close()

    alignment = [seq for seq in SeqIO.parse(proc.stdout, 'fasta')]

    with open(out_fasta_path, 'w') as f:
        SeqIO.write(alignment, f, 'fasta')

    #Check if there is any kdr mutation listed in kdr_list.json

    kdr_list = os.path.dirname(os.path.abspath(__file__)) + "/../scripts/kdr_list.json"
    with open(kdr_list) as f:
        kdr_dict = json.load(f)

    each_pos = [0, 0, 0]
    for pos in range(len(alignment[0].seq)):
        for seq in [0, 1, 2]:
            if alignment[seq][pos] != "-":
                each_pos[seq] += 1
        AA_Change_ck = ""
        AA_Change_dl = ""
        is_kdr_ck = ""
        is_kdr_dl = ""
        if alignment[0][pos] != alignment[1][pos] and\
                    not (alignment[0][pos] == "-" or alignment[1][pos] == "-"):
            AA_Change_ck = "{}{}".format(each_pos[0], alignment[1][pos])
            if AA_Change_ck in kdr_dict:
                if kdr_dict[AA_Change_ck] == 1:
                    is_kdr_ck = "kdr!"
                else:
                    is_kdr_ck = "kdr?"
        if alignment[0][pos] != alignment[2][pos] and\
                    not (alignment[0][pos] == "-" or alignment[2][pos] == "-"):
            AA_Change_dl = "{}{}".format(each_pos[0], alignment[2][pos])
            if AA_Change_dl in kdr_dict:
                if kdr_dict[AA_Change_dl] == 1:
                    is_kdr_dl = "kdr!"
                else:
                    is_kdr_dl = "kdr?"

        if AA_Change_ck and AA_Change_ck == AA_Change_dl:
            print("    {} in ck and dl {}".format(AA_Change_ck, is_kdr_ck), file=info)
        elif AA_Change_ck:
            print("    {} in ck {}".format(AA_Change_ck, is_kdr_ck), file=info)
        elif AA_Change_dl:
            print("    {} in dl {}".format(AA_Change_dl, is_kdr_dl), file=info)

    info.close()

if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description=
                            "Make AA aligment with M. domestica VGSC"
                            "from dna fasta and bed")

    parser.add_argument("ref_fa",
                       help = 'Reference fasta')

    parser.add_argument("bed",
                       help = "Reference annotation bed")

    parser.add_argument("-o", "--out_aln_fasta", dest = "out_fasta_path",
                       help = "Multiple AA aligment fasta file to output.")

    parser.add_argument("-t", "--translate", dest = "translate",
                       help = "Output only translation to this file.")

    parser.add_argument("-m", "--mdom_path", dest = "mdom_path",
                       help = 'M. domestica aa fasta path [MoNaS/misc/AAB47604.fa]',
                       default = script_dir + "/../misc/AAB47604.fa")

    args = parser.parse_args()

    if not (args.out_fasta_path or args.translate):
        sys.exit(print("Error!: Neither -o nor -t has been given."))

    for file in [args.ref_fa, args.bed, args.mdom_path]:
        if not os.path.isfile(file):
            print(file + " was not found!", file = sys.stderr)
            sys.exit(1)

    main(ref_path = args.ref_fa,
         bed_path = args.bed,
         Mdom_path = args.mdom_path,
         out_fasta_path = args.out_fasta_path,
         translate = args.translate)

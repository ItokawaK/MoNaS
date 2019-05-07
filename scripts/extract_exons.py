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
import os
from finalize_table import Bed

def extract(fa_path, bed_path, flanking = 0):

    bed = Bed(bed_path)

    out = []

    for seq in SeqIO.parse(fa_path, "fasta"):
        if seq.id == bed.chrom:
            if bed.strand == "+":
                bed.entries.sort(key = lambda ex: ex.start)
            else:
                bed.entries.sort(key = lambda ex: ex.start, reverse = True)
            for bed_l in bed.entries:
                exon = seq[(bed_l.start - flanking):(bed_l.end + flanking)]
                if bed.strand == '+':
                    exon.id = bed_l.name
                    out.append(exon)
                else:
                    exon = exon.reverse_complement()
                    exon.id = bed_l.name
                    out.append(exon)

    return out

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description='Genotype VGSC gene.')

    parser.add_argument("ref_fasta",
                       help = 'Reference genome fasta')

    parser.add_argument("ref_bed",
                       help = "Reference BED file")

    parser.add_argument("-f", "--flanking",
                       help = "Extract this additional bases of flannking ends",
                       type = int,
                       default = 0)

    # parser.add_argument("-o", "--out_fasta", dest = "out_fasta_path",
    #                    help = "fasta to output")

    args = parser.parse_args()

    fa_path = args.ref_fasta
    bed_path = args.ref_bed
    flanking = args.flanking

    for file in [fa_path, bed_path]:
        if not os.path.isfile(file):
            print(file + " was not found!", file = sys.stderr)
            sys.exit(1)

    for seq in extract(fa_path, bed_path, flanking):
        print(">" + seq.id)
        print(str(seq.seq))

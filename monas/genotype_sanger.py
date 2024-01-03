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

import os
import sys
import subprocess
import math
from Bio import SeqIO
import shutil

def chop_dna(dna):
    """ Chopping dna into small (150 bp) peices"""
    read_len = 150
    max_ovl = 50
    min_coverage = 5
    out = []

    dna_len = len(dna)
    base_id = dna.id
    starts = []
    start = 0
    read_n = math.floor((dna_len - max_ovl)/(read_len - max_ovl))
    if read_n > 1:
        ovl_len = (read_len * read_n - dna_len)/(read_n - 1)
    else:
        ovl_len = max_ovl

    cnt = 0
    for i in range(read_n):
        for ii in range(min_coverage):
            if i == read_n - 1:
                out_seq = dna[int(start) : ]
            else:
                out_seq = dna[int(start) : int(start + read_len)]

            out_seq.id = base_id + "_" + str(cnt)
            out_seq.letter_annotations["phred_quality"] = [40] * len(out_seq)
            out.append(out_seq)
            cnt += 1

        start += (read_len - ovl_len)

    return out

def write_fasta(in_seq_fasta, out_fa_dir, out_list):
    """Out put fasta file from given seq list"""

    if not os.path.isdir(out_fa_dir):
        os.mkdir(out_fa_dir)

    seqs = SeqIO.parse(in_seq_fasta, 'fasta')
    out_list_lines = []
    for seq in seqs:
        chopped_dna = chop_dna(seq)
        out_file_path = os.path.join(out_fa_dir, seq.id + ".fq")
        out_list_lines.append(seq.id + " " + os.path.abspath(out_file_path))
        with open(out_file_path, 'w') as f:
            SeqIO.write(chopped_dna, f, 'fastq')

    with open(out_list, 'w') as f:
        for l in out_list_lines:
            f.write(l + '\n')

def main():
    import argparse
    import random
    import atexit

    parser = argparse.ArgumentParser(description = 'Find mutations on VGSC from sangar sequences')

    parser.add_argument('-s', '--species', dest = 'species',
                        help = 'Species name. It should be same to the dirname of references.')
    parser.add_argument('-r', '--ref_root', dest = 'ref_root',
                        help = 'Root directly of references. Deault = MoNaS/references')
    parser.add_argument('-t', '--max_cpu', dest = 'num_cpu',
                       type = int,
                       default = 4,
                       help = 'Maximum number of threads. [4]')
    parser.add_argument('-b', '--bwa_treads', dest = 'num_threads',
                       type = int,
                       default = 1,
                       help = 'Number of treads per bwa process. [1]')
    parser.add_argument('-o', '--out_table', dest = 'out_tbl',
                       help = 'Name of the out table. Should be new. [out_table.tsv]',
                       default  = "out_table.tsv")
    parser.add_argument('-n', '--no_clean', dest = 'do_clean',
                         action ='store_false',
                         default = True,
                         help = 'Do not clean temp files. [off]'
                         )
    parser.add_argument('in_fasta',
                         help = 'Input sequences in multi fasta format.')

    args = parser.parse_args()

    if not os.path.isfile(args.in_fasta):
        print("Error!: " + args.in_fasta + " does not exist.",
              file = sys.stderr)
        sys.exit(1)
    else:
        in_seq_fasta = args.in_fasta

    if not args.species:
        print("Error!: -s option is mandately.", file = sys.stderr)
        sys.exit(1)
    else:
        species = args.species

    if os.path.isfile(args.out_tbl):
        print("Error!: " + args.out_tbl + " already exists.",
              file = sys.stderr)
        sys.exit(1)
    else:
        out_tbl = args.out_tbl

    rndm = "".join(random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                  "abcdefghijklmnopqrstuvwxyz"
                                  "0123456789", k=15))
    out_fa_dir = rndm + "_fastq_tmp"
    out_list = rndm + ".list.tmp"
    monas_out_path = rndm + "_monasout_tmp"

    def cleaner(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)
        else:
            print(path + " was not found.", file=sys.stderr)

    if args.do_clean:
        atexit.register(cleaner, monas_out_path)
        atexit.register(cleaner, out_list)
        atexit.register(cleaner, out_fa_dir)

    write_fasta(in_seq_fasta, out_fa_dir, out_list)

    monas_path = os.path.dirname(os.path.abspath(__file__)) + "/genotype.py"

    cmd = [monas_path,
          "-s", species,
          "-l", out_list,
          "-o", monas_out_path,
          "-t", str(args.num_cpu),
          "-b", str(args.num_threads)]

    if args.ref_root:
        cmd += ["-r", args.ref_root]

    subprocess.call(cmd)

    if os.path.isfile(monas_out_path + "/table_with_Mdomcoord.tsv"):
        shutil.copy(monas_out_path + "/table_with_Mdomcoord.tsv", out_tbl)
    else:
        print("Error!", file=sys.stderr)

    # if args.do_clean:
    #     shutil.rmtree(monas_out_path)
    #     shutil.rmtree(out_fa_dir)
    #     os.remove(out_list)

if __name__ == '__main__':
    main()
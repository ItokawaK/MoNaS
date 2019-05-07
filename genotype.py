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
import shutil
import subprocess
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from scripts import finalize_table
from scripts.configuration import GenomeRef
from scripts.jobs import Job

version = "1.0"

def parse_sample_list(sample_list):
    # Read sample list file
    # Returns a list of lists [sample_name, fq.gz path1, fq.gz path2]

    samples = []
    with open(sample_list) as f:
        line = f.readline()
        while line:
            line = line.strip()
            samples.append(line.split(" "))
            line = f.readline()

    return(samples)

def usage():
    return (
    """
                  __  __       _   _        _____
                 |..\\/..|     |.\\ |.|      /.____|
                 |.\\../.| ___ |..\\|.| __ _|.(___
                 |.|\\/|.|/._.\\|...`.|/._`.|\\___ \\
                 |.|  |.|.(_).|.|\\..|.(_|.|____).|
                 |_|  |_|\\___/|_| \\_|\\__,_|_____/

    genotyp.py -s species_name -o out_dir_path -l list_file_path
                 (-t num_max_threads -b num_threads_per_bwa
                  -m mode[ngs_dna|ngs_rna] -r ref_dir_path
                  -v variant_caller[freebayes|gatk]
                  )
    """)

def description(version):
    return (
       " MoNaS (version {}) - A program genotyping VGSC genes from NGS reads.".format(version)
     )
if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description = description(version),
                                     usage = usage())

    parser.add_argument('-s', '--species', dest = 'species',
                        help = 'Species name. It should be same to the dirname of references.')
    parser.add_argument('-l', '--sample_list', dest = 'sample_list',
                       help = 'Path for a list file decribing name of sampeles and fastq files.')
    parser.add_argument('-t', '--max_cpu', dest = 'num_cpu',
                       type = int,
                       default = 4,
                       help = 'Maximum number of threads. [4]')
    parser.add_argument('-b', '--bwa_treads', dest = 'num_threads',
                       type = int,
                       default = 4,
                       help = 'Number of treads per bwa process. [4]')
    parser.add_argument('-o', '--out_dir', dest = 'out_dir',
                       help = 'Name of out directly. Should be new.')
    parser.add_argument('-r', '--ref_root', dest = 'ref_root',
                        default = script_dir + "/references",
                        help = 'Root directly of references. Deault = MoNaS/references')
    parser.add_argument('-m', '--mode', dest = 'mode',
                         default = 'ngs_dna',
                         choices = ['ngs_dna', 'ngs_rna'],
                         help = 'Analysis mode. [ngs_dna]'
                         )
    parser.add_argument('-c', '--variant_caller', dest = 'variant_caller',
                         default = "freebayes",
                         choices = ['freebayes', 'gatk'],
                         help = 'Variant caller to be used. Default is freebayes.'
                         )
    parser.add_argument('-n', '--no_clean', dest = 'do_clean',
                         action='store_false',
                         default = True,
                         help = 'Do not clean old BAM files after rmdup. Off in default.'
                         )
    parser.add_argument('-v', '--version', dest = 'show_version',
                        action='store_true',
                        default = False,
                        help = 'Show version and exit.'
                        )

    args = parser.parse_args()

    if args.show_version:
        print(version)
        sys.exit(0)

    if not (args.species and args.out_dir and args.sample_list):
        print("Error: Species, out_dir and list are mandately!", file = sys.stderr)
        print("USAGE", file = sys.stderr)
        print(usage(), file = sys.stderr)
        sys.exit(1)

    num_cpu = args.num_cpu
    if args.num_threads:
        num_threads = args.num_threads
        if  num_threads < num_cpu:
            num_proc = num_cpu // num_threads
        else:
            num_threads = num_cpu
            num_proc = 1
    else:
        if num_cpu > 3:
            num_threads = 4
        else:
            num_threads = 1
        num_proc = num_cpu // num_threads

    out_dir = args.out_dir
    out_bam_dir1 = out_dir + "/BAMs"
    out_bam_dir2 = out_dir + "/BAMs_rmdup"
    vcf_out_dir = out_dir + "/VCFs"
    #bin_path = Bin(args.bin_root) #Bin object
    genome_ref = GenomeRef(args.ref_root, args.species)  #GenomeRef object

    #configuration for genotyping
    genome_ref.check_program_path(args.mode, args.variant_caller)
    genome_ref.check_genomedb(args.mode, args.variant_caller, num_cpu = args.num_cpu)
    genome_ref.check_existence_for_genotype()


    out_table = out_dir + "/out_table"

    if os.path.isfile(args.sample_list):
        sample_list = args.sample_list
        samples = parse_sample_list(sample_list)
    else:
        print(args.sample_list + r" does not exist /(*o*)\ ")
        sys.exit(1)

    if os.path.exists(out_dir):
        print("Error: " + out_dir + " already exists!", file = sys.stderr)
        sys.exit(1)

    job = Job(genome_ref, args.mode)


    #Create output dir
    os.mkdir(out_dir)

    if not os.path.isdir(out_bam_dir1):
        os.mkdir(out_bam_dir1)

    job.map_and_sort_mp(num_threads = num_threads,
                        num_proc = num_proc,
                        samples = samples,
                        out_bam_dir = out_bam_dir1)

    if not os.path.isdir(out_bam_dir2):
        os.mkdir(out_bam_dir2)

    job.rmdup_and_index_mp(num_cpu = num_cpu,
                           in_bams = job.bams_to_process,
                           out_bam_dir = out_bam_dir2)

    if args.do_clean:
        shutil.rmtree(out_bam_dir1)

    if args.variant_caller == "gatk":
        if not os.path.isdir(vcf_out_dir):
            os.mkdir(vcf_out_dir)
        csqvcfs = job.variant_analysis_gatk_mp(num_cpu = 12,
                                               in_bams = job.bams_to_process,
                                               vcf_out_dir = vcf_out_dir)
    else:
        job.variant_analysis_fb(num_cpu = num_cpu,
                                in_bams = job.bams_to_process,
                                out_vcf = out_dir + "/out.vcf",
                                out_csqvcf = out_dir + "/out_csq.vcf")

        csqvcfs = [out_dir + "/out_csq.vcf"]

    finalize_table.create_table(
                     csqvcfs = csqvcfs,
                     info_to_get = ["CHROM",
                                    "POS",
                                    "REF_ALLELE",
                                    "ALT_ALLELE",
                                    "QUAL",
                                    "GENOTYPE",
                                    "AA_CHANGE",
                                    "AA_CHANGE_MDOM",
                                    "KDR_EVIDENCE",
                                    "AD",
                                    "EXON"],
                     # header = "#ID\tCHROM\tPOS\tREF_ALLELE\tALT_ALLELE\tQUAL\tGT\t"
                     #          "AA_CHANGE\tAA_CHANGE_HOUSEFLY\tAD\tEXON",
                     bed_file = genome_ref.bed,
                     mdom_fasta = genome_ref.mdom_fa,
                     out_table_file = out_dir + "/table_with_Mdomcoord.tsv"
                     )

    print(r"MoNaS is done \(^o^)/ !", file = sys.stderr)

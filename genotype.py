#!/usr/bin/env python
# coding: utf-8

import os
import shutil
import subprocess
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from scripts import finalize_table
from scripts.configuration import GenomeRef
from scripts.jobs import Job

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
    return ("genotyp.py -s species -o out_dir\n"
            "\t         (--max_cpu int --ref_root path\n"
            "\t          --sample_list path --fasta path\n"
            "\t          --mode ngs_dna|ngs_rna|sanger_dna)\n")

if __name__ == '__main__':

    script_dir = sys.path[0]

    parser = argparse.ArgumentParser(description = 'Genotype VGSC gene.',
                                     usage = usage())

    parser.add_argument('-s', '--species', dest = 'species',
                        help = 'Species name. It should be same to the dirname of references.')
    parser.add_argument('-l', '--sample_list', dest = 'sample_list',
                       help = 'Path for a list file decribing name of sampeles and fastq files.')
    parser.add_argument('-t', '--max_cpu', dest = 'num_cpu',
                       type = int,
                       default = 4,
                       help = 'Maximum number of threads. [4]')
    parser.add_argument('-o', '--out_dir', dest = 'out_dir',
                       help = 'Name of out directly. Should be new.')
    parser.add_argument('-r', '--ref_root', dest = 'ref_root',
                        default = script_dir + "/references",
                        help = 'Root directly of references. Deault = MoNaS_v1.x/references')
    parser.add_argument('-f', '--fasta', dest = 'fasta',
                         help = 'Path to fasta file of Sanger seq.')
    parser.add_argument('-m', '--mode', dest = 'mode',
                         default = 'ngs_dna',
                         choices = ['ngs_dna', 'ngs_rna', 'sanger_dna'],
                         help = 'Analysis mode. [ngs_dna]'
                         )
    parser.add_argument('-v', '--variant_caller', dest = 'variant_caller',
                         default = "freebayes",
                         choices = ['freebayes', 'gatk'],
                         help = 'Variant caller to be used. Default is freebayes.'
                         )
    parser.add_argument('-n', '--no_clean', dest = 'do_clean',
                         action='store_false',
                         default = True,
                         help = 'Do not clean old BAM files after rmdup. Off in default.'
                         )



    args = parser.parse_args()

    if not (args.species and args.ref_root and args.out_dir):
        print("Error: species and out_dir are mandately!", file = sys.stderr)
        print(usage())
        sys.exit(1)

    if args.mode in ['ngs_dna', 'ngs_rna']:
        if not args.sample_list:
            print("Error: sample_list is mandately for NGS data!", file = sys.stderr)
            print(usage())
            sys.exit(1)
    else:
        if not args.fasta:
            print("Error: fasta is mandately for Sanger data!", file = sys.stderr)
            print(usage())
            sys.exit(1)

    num_cpu = args.num_cpu
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
    genome_ref = GenomeRef(args.ref_root,
                           args.species,
                           args.mode,
                           args.num_cpu,
                           args.variant_caller) #GenomeRef object
    out_table = out_dir + "/out_table"

    if os.path.isfile(args.sample_list):
        sample_list = args.sample_list
        samples = parse_sample_list(sample_list)
    else:
        print(args.sample_list + r" does not exist /(*o*)\ ")
        sys.exit(1)

    if os.path.exists(out_dir):
        print("Warning: " + out_dir + " already exists!", file = sys.stderr)
        print("  Do you wanna proceed anyway (._.)? 0: no, >0: yes", file = sys.stderr)
        to_proceed = input(">>>>  ")
        if to_proceed == "0":
            sys.exit(1)

    job = Job(genome_ref)


    #Create output dir
    if not os.path.isdir(out_dir):
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
        job.variant_analysis_gatk_mp(num_cpu = 12,
                                     in_bams = job.bams_to_process,
                                     vcf_out_dir = vcf_out_dir,
                                     out_table = out_table)
    else:
        job.variant_analysis_fb(in_bam_list = job.bams_to_process,
                                out_vcf = out_dir + "/out.vcf",
                                out_csqvcf = out_dir + "/out_csq.vcf",
                                out_table = out_table)

    finalize_table.create_table(out_dir + "/out_csq.vcf",
                                genome_ref.bed,
                                genome_ref.mdom_fa,
                                script_dir + "/scripts/kdr_list.json",
                                out_dir + "/table_with_Mdomcoord.tsv")

    print(r"MoNaS is done \(^o^)/ !", file = sys.stderr)

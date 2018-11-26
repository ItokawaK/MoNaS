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

    species = args.species
    num_cpu = args.num_cpu
    if num_cpu > 3:
        num_threads = 4
    else:
        num_threads = 1
    num_proc = num_cpu // num_threads
    out_dir = args.out_dir
    out_bam_dir = out_dir + "/BAMs"
    vcf_out_dir = out_dir + "/VCFs"
    #bin_path = Bin(args.bin_root) #Bin object
    ref_dir = args.ref_root
    genome_ref = GenomeRef(ref_dir, species) #GenomeRef object
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
    if not os.path.isdir(out_bam_dir):
        os.mkdir(out_bam_dir)

    #Executing multiprocesses of bwa mem | samtools sort
    with ProcessPoolExecutor(max_workers = num_proc) as executor:
        executed = [executor.submit(job.map_and_sort,
                                    sample,
                                    num_threads,
                                    out_bam_dir) for sample in samples]

    out_bams = [ex.result() for ex in executed]

    #Executing multiprocesses of samtools rmdup and samtools index
    with ProcessPoolExecutor(max_workers = num_cpu) as executor:
        executed = [executor.submit(job.rmdup_and_index, bam) for bam in out_bams]

    if not os.path.isdir(vcf_out_dir):
        os.mkdir(vcf_out_dir)

    #Executing multiprocesses of gatk
    #Output files include table queried by bcftools
    executed2 = []
    with ProcessPoolExecutor(max_workers = 12) as executor:
        for bam in out_bams:
            vcf_name = os.path.basename(bam).rstrip(".bam")
            executed2.append(executor.submit(job.variant_analysis,
                                             vcf_out_dir,
                                             bam,
                                             vcf_name)
                            )

    single_tables = [ex.result() for ex in executed2]

    #Concatenating tables
    with open(out_table, 'w') as outfile:
        for single_table in single_tables:
            with open(single_table) as infile:
                outfile.write(infile.read())

    finalize_table.create_table(out_table,
                                genome_ref.bed,
                                genome_ref.mdom_fa,
                                out_dir + "/table_with_Mdomcoord.tsv")

    print(r"MoNaS is done \(^o^)/ !", file = sys.stderr)

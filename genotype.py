#!/usr/bin/env python
# coding: utf-8

import os
import shutil
import subprocess
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
import finalize_table
from configuration import GenomeRef
from jobs import Job

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

if __name__ == '__main__':

    script_dir = sys.path[0]

    parser = argparse.ArgumentParser(description = 'Genotype VGSC gene.')

    parser.add_argument("species",
                       choices = ['Cpip', 'Aalb', 'Aaeg'],
                       help = 'Species name. It should be same to the dirname of references.')
    parser.add_argument("sample_list",
                       help = 'Path for a list file decribing name of sampeles and fastq files.')
    parser.add_argument("num_cpu",
                       type = int,
                       help = 'Number of threads to use.')
    parser.add_argument("out_dir",
                       help = 'Name of out directly. Should be new.')
    parser.add_argument('-r', '--ref_root', dest = 'ref_root',
                        default = script_dir + "/references",
                        help = 'Root directly of references. Deault = MoNaS_v1.x/references')
    parser.add_argument('-b', '--bin_root', dest = 'bin_root',
                        default = script_dir + "/bin",
                        help = 'Root directly of third party binaries. Deault = MoNaS_v1.x/bin')


    args = parser.parse_args()

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
        print("Warning: " + out_dir + " already exists!")
        print("  Do you wanna proceed anyway (._.)? 0: no, >0: yes")
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

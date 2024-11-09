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
import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
import re

try:
    from monas import finalize_table
    from monas.configuration import GenomeRef
    from monas.jobs import Job
    from monas import logging_conf
    from monas.__version__ import __version__
except:
    import finalize_table
    from configuration import GenomeRef
    from jobs import Job
    import logging_conf
    from __version__ import __version__

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

def check_versions():

    versions ={}    
    versions['bwa'] = 'Unknown'
    versions['samtools'] = 'Unknown'
    versions['freebayes'] = 'Unknown'

    # bwa
    sbp = subprocess.run('bwa', stderr=subprocess.PIPE, universal_newlines=True)
    for l in sbp.stderr.split('\n'):
        match = re.findall('Version: (.+)', l)
        if match:
            versions['bwa'] = match[0]
            break

    # samtools 
    sbp = subprocess.run(['samtools', 'version'], stdout=subprocess.PIPE, universal_newlines=True)
    versions['samtools'] = sbp.stdout.split('\n')[0].split(' ')[1]

    # freeBayes
    sbp = subprocess.run(['freebayes'], stdout=subprocess.PIPE, universal_newlines=True)
    for l in sbp.stdout.split('\n'):
        match = re.findall('version: +(.+)', l)
        if match:
            versions['freebayes'] =  match[0]
            break

    return versions

def usage():
    return (
    """
                  __  __       _   _        _____
                 |..\\/..|     |.\\ |.|      /.____|
                 |.\\../.| ___ |..\\|.| __ _|.(___
                 |.|\\/|.|/._.\\|...`.|/._`.|\\___ \\
                 |.|  |.|.(_).|.|\\..|.(_|.|____).|
                 |_|  |_|\\___/|_| \\_|\\__,_|_____/

    genotype.py -s species_name -o out_dir_path -l list_file_path
                 (-t num_max_threads -b num_threads_per_bwa
                  -m mode[ngs_dna|ngs_rna] -r ref_dir_path
                  -c variant_caller[freebayes|gatk],
                  --version
                  )
    """)

# def description(version):
#     return (
#        " MoNaS (version {}) - A program genotyping VGSC genes from NGS reads.".format(version)
#      )

def main(args):

    if not (args.species and args.out_dir and args.sample_list):
        print("Error: Species, out_dir and list are mandatory!", file = sys.stderr)
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
    out_table =  out_dir + "/table_with_Mdomcoord.tsv"
    out_stats =  out_dir + "/stats"


    if os.path.isfile(args.sample_list):
        sample_list = args.sample_list
        samples = parse_sample_list(sample_list)
    else:
        print(args.sample_list + r" does not exist /(*o*)\ ")
        sys.exit(1)

    if(len(samples) == 0):
        print("Error: There was no sample.")
        sys.exit(1)
    else:
        file_not_found = []
        for sample in samples:
            tmp = []
            for i in (1,2):
                if not os.path.isfile(sample[i]):
                    file_not_found.append(sample[i])
        if(len(file_not_found) > 0):
            print("Error: Follwing files were not found.")
            print(file_not_found)
            sys.exit(1)

    #Create output dir
    if os.path.exists(out_dir) and not args.resume:
        print("Error: " + out_dir + " already exists!", file = sys.stderr)
        sys.exit(1)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    #----logging setting----
    log_file = os.path.join(args.out_dir, 'log')
    logger = getLogger(__name__)
    logging_conf.set(log_file)
    logger.info("Aanalysis started")
    logger.info(f'MoNaS version: {__version__}')
    other_app_versions = check_versions()
    for k,v in other_app_versions.items():
        logger.info(f'{k} version: {v}')
    #----end----

    #bin_path = Bin(args.bin_root) #Bin object
    genome_ref = GenomeRef(args.ref_root, args.species)  #GenomeRef object

    #configuration for genotyping
    genome_ref.check_program_path(args.mode, args.variant_caller)
    genome_ref.check_genomedb(args.mode, args.variant_caller, num_cpu = args.num_cpu)
    genome_ref.check_existence_for_genotype()

    job = Job(genome_ref, args.mode, out_dir)

    if not os.path.isdir(out_bam_dir1):
        os.mkdir(out_bam_dir1)

    logger.info("[1/2] Start mapping and sorting for {} samples...".format(len(samples)))
    logger.info("  {} BWA processes are running concurrently, "
                "{} therads are used at each process".format(num_proc, num_threads)
                )

    job.map_and_sort_mp(num_threads = num_threads,
                        num_proc = num_proc,
                        samples = samples,
                        out_bam_dir = out_bam_dir1)

    logger.info("...Finished mapping and sorting for all samples.")
    logger.info("Calculating mapping stats...")

    job.get_stats_mp(num_threads=num_threads,
                     samples = samples,
                     in_bam_dir = out_bam_dir1,
                     out_file = out_stats
                     )

    # if not os.path.isdir(out_bam_dir2):
    #     os.mkdir(out_bam_dir2)
    #
    # logger.info("[2/3] Start removing PCR duplicates...")

    # job.rmdup_and_index_mp(num_cpu = num_cpu,
    #                        in_bams = job.bams_to_process,
    #                        out_bam_dir = out_bam_dir2)
    #
    # logger.info("...Finished removing PCR duplicates.")
    #
    # if args.do_clean:
    #     shutil.rmtree(out_bam_dir1)

    logger.info("[2/2] Started variant calling using {}... ".format(args.variant_caller))

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

    logger.info("...Finished variant calling")

    logger.info("Writing the final result in table_with_Mdomcoord.tsv...")
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
                     out_table_file = out_table
                     )

    logger.info(r"MoNaS is finished!")

# if __name__ == '__main__':
#     main()

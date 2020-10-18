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
import subprocess
import sys
import shutil
import logging

from scripts.configuration import GenomeRef
from concurrent.futures import ProcessPoolExecutor

class Job:
    def __init__(self, genome_ref, mode, out_dir):
        self.gref = genome_ref
        self.bams_to_process = [] #List of bam files to be processed in next step
        self.mode = mode
        self.logger = logging.getLogger(__name__)
        self.out_dir = out_dir
        self.log_file = os.path.join(self.out_dir, "job.err")

    def map_and_sort(self, sample, num_threads, out_bam_dir):
        #Mapping and sort
        # 'sample' sould be a list of [name, fastq1_path, fastq2_path]
        # retruns paths of resulted bam files in 'out_bam_dir'

        out_bam_path = out_bam_dir + "/" + sample[0] + ".bam"

        if self.mode in ['ngs_dna']:
            cmd1 = ['bwa', "mem",
                    "-t", str(num_threads),
                    "-R", r'@RG\tID:ID_' + sample[0] + r'\tSM:' + sample[0],
                    self.gref.bwa_db] + sample[1:3]
        elif self.mode in ['ngs_rna']:
            if sample[2]:
                sample_cmd = ["-1", sample[1], "-2", sample[2]]
            else:
                sample_cmd = ["-r", sample[1]]

            cmd1 = ["hisat2",
                    "-x", self.gref.hisat_db,
                    "-p", str(num_threads),
                    "--rg-id", r'ID_' + sample[0],
                    "--rg", "SM:" + sample[0],
                    ] + sample_cmd

        cmd2 = ['samtools', "sort",
                "-o", out_bam_path]

        # Uncomment to skip if bam.bai already exists. For debuggin variant calling.
        #if os.path.isfile(out_bam_path + ".bai"):
        #    return(out_bam_path)

        with open(self.log_file, "a") as err_log:
            try:
                proc1 = subprocess.Popen(
                             cmd1,
                             stdout = subprocess.PIPE,
                             stderr = err_log)
                proc2 = subprocess.Popen(
                             cmd2,
                             stdin = proc1.stdout,
                             stderr = err_log)
                proc2.communicate()

                self.logger.info("   Finished mapping and sorting: " + sample[0])
                return(out_bam_path)
            except:
                self.logger.info("   Failed mapping and sorting: " + sample[0])

    def map_and_sort_mp(self, num_threads, num_proc, samples, out_bam_dir):
        with  ProcessPoolExecutor(max_workers = num_proc) as executor:
            executed = [executor.submit(self.map_and_sort,
                                        sample,
                                        num_threads,
                                        out_bam_dir) for sample in samples]

        self.bams_to_process = [ex.result() for ex in executed]

    def rmdup_and_index(self, in_bam_path, out_bam_dir):
        # samptools rmdup and index for a given bam
        # Result is stored in same dir
        # Original bam file will be relaced

        in_bam_dir = os.path.dirname(in_bam_path)
        in_bam_basename = os.path.basename(in_bam_path)
        out_bam_path = os.path.join(out_bam_dir, in_bam_basename)

        cmd1 = ['samtools', "rmdup",
               in_bam_path,
               out_bam_path
               ]

        cmd2 = ['samtools', "index",
               out_bam_path
               ]

        with open(self.log_file, "a") as err_log:
            try:
                p = subprocess.Popen(
                     cmd1,
                     stdout = err_log,
                     stderr = err_log).communicate()
                p = subprocess.Popen(
                     cmd2,
                     stdout = err_log,
                     stderr = err_log).communicate()
                return(out_bam_path)
            except:
                self.logger.error("   Failed rmdup and indexing: " + sample[0])

    def rmdup_and_index_mp(self, num_cpu, in_bams, out_bam_dir):
        with ProcessPoolExecutor(max_workers = num_cpu) as executor:
            executed = [executor.submit(self.rmdup_and_index,
                                        bam,
                                        out_bam_dir
                                        ) for bam in in_bams]

        self.bams_to_process = [ex.result() for ex in executed]

    def variant_analysis_gatk(self, out_dir, in_bam, sample_name):
        # Variant calling using gatk and annotion usig bcftools csq
        # Returns path of out_table

        out_gvcf = out_dir + "/" + sample_name + ".g.vcf"
        out_vcf = out_dir + "/" + sample_name + ".vcf"
        out_csqvcf = out_dir + "/" + sample_name + "_csq.vcf"
        out_table = out_dir + "/" + sample_name + "_out_table"

        cmd1 = ['gatk', "HaplotypeCaller",
               "-L", self.gref.bed,
               "-R", self.gref.ref_fa,
               "-I", in_bam,
               "-O", out_gvcf,
               "-ERC", "GVCF",
               "--max-mnp-distance", "5",
               "--native-pair-hmm-threads", "1"]

        cmd2 = ['gatk', "GenotypeGVCFs",
               "-R", self.gref.ref_fa,
               "-V", out_gvcf,
               "-O", out_vcf]

        cmd3 = ['bcftools', "csq",
               "-g", self.gref.gff3,
               "-f", self.gref.ref_fa,
               "-p", "a",
               "-l",
               "-o", out_csqvcf,
               out_vcf]

        with open(self.log_file, "a") as err_log:
            p = subprocess.Popen(cmd1, stderr=err_log).communicate()
            p = subprocess.Popen(cmd2, stderr=err_log).communicate()
            p = subprocess.Popen(cmd3, stderr=err_log).communicate()

        return(out_csqvcf)

    def variant_analysis_gatk_mp(self, num_cpu, in_bams, vcf_out_dir):

        executed2 = []

        with ProcessPoolExecutor(max_workers = num_cpu) as executor:
            for bam in in_bams:
                vcf_name = os.path.basename(bam).rstrip(".bam")
                executed2.append(executor.submit(self.variant_analysis_gatk,
                                                 vcf_out_dir,
                                                 bam,
                                                 vcf_name)
                                )

        csqvcfs = [ex.result() for ex in executed2]

        return csqvcfs


    def run_freebayes(self, in_bams, region):
        cmd = ["freebayes",
                   "-r", region,
                   "-f", self.gref.ref_fa] + in_bams

        with open(self.log_file, "a") as err_log:
            proc = subprocess.Popen(cmd,
                      stdout=subprocess.PIPE,
                      stderr=err_log)

        header = []
        body = []

        stdout_data = proc.communicate()[0]
        lines = stdout_data.decode().split("\n")
        for line in lines:
            line = line
            if line.startswith("#"):
                header.append(line)
            else:
                body.append(line)

        return([header, body])


    def variant_analysis_fb(self, num_cpu, in_bams, out_vcf, out_csqvcf):

        bed = []
        with open(self.gref.bed) as f:
            bed = [l.rstrip().split("\t") for l in f.readlines()]

        regions = [b[0] + ":" + b[1] + "-" + b[2] for b in bed if b != '']

        self.logger.info("   Conducting for {} partioned regions...".format(len(regions)))
        with ProcessPoolExecutor(max_workers = num_cpu) as executor:
            executed = [executor.submit(self.run_freebayes,
                                        in_bams,
                                        region) for region in regions]

        header = executed[0].result()[0]
        bodies = [ex.result()[1] for ex in executed]

        with open(out_vcf, 'w') as f:

            for l in header:
                if l.rstrip() != '':
                   print(l.rstrip(), file = f)
            for body in bodies:
                for l in body:
                    if l.rstrip() != '':
                        print(l.rstrip(), file = f)

        cmd = ['bcftools', "csq",
                "-g", self.gref.gff3,
                "-f", self.gref.ref_fa,
                "-p", "a",
                "-l",
                "-o", out_csqvcf,
                out_vcf]

        with open(self.log_file, "a") as err_log:
            p = subprocess.Popen(cmd, stderr=err_log).communicate()

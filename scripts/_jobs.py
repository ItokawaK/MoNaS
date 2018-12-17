import os
import subprocess
import sys
from scripts.configuration import GenomeRef
from concurrent.futures import ProcessPoolExecutor
import shutil

class Job:
    def __init__(self, genome_ref):
        self.gref = genome_ref

    def map_and_sort(self, sample, num_threads, out_bam_dir):
        #Mapping and sort
        # 'sample' sould be a list of [name, fastq1_path, fastq2_path]
        # retruns paths of resulted bam files in 'out_bam_dir'

        out_bam_path = out_bam_dir + "/" + sample[0] + ".bam"

        if self.gref.mode in ['ngs_dna']:
            cmd1 = ['bwa', "mem",
                    "-t", str(num_threads),
                    "-R", r'@RG\tID:ID_' + sample[0] + r'\tSM:' + sample[0],
                    self.gref.bwa_db] + sample[1:3]
        elif self.gref.mode in ['ngs_rna']:
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

        proc1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        proc2 = subprocess.Popen(cmd2, stdin = proc1.stdout)
        proc2.communicate()
        return(out_bam_path)

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

        subprocess.call(cmd1)
        subprocess.call(cmd2)

        return(out_bam_path)

    def variant_analysis(self, out_dir, in_bam, sample_name):
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
               "-p", "a"
               "-l",
               "-o", out_csqvcf,
               out_vcf]

        proc1 = subprocess.call(cmd1)
        proc2 = subprocess.call(cmd2)
        proc3 = subprocess.call(cmd3)

        cmd4 = ['bcftools', "query",
               "-f", r"[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TBCSQ\t%TGT\t%AD\n]",
               "-o", out_table,
               out_csqvcf]

        subprocess.call(cmd4)
        return(out_table)

        #with open(final_out_table, 'w') as f:
        #    subprocess.Popen(cmd4, stdout = f)
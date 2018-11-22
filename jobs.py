import os
import subprocess
import sys
from configuration import GenomeRef
from configuration import Bin
from concurrent.futures import ProcessPoolExecutor

class Job:
    def __init__(self, bin, genome_ref):
        self.bin = bin # Bin object
        self.gref = genome_ref

    def map_and_sort(self, sample, num_threads, out_bam_dir):
        #Mapping and sort
        # 'sample' sould be a list of [name, fastq1_path, fastq2_path]
        # retruns paths of resulted bam files in 'out_bam_dir'

        out_bam_path = out_bam_dir + "/" + sample[0] + ".bam"

        cmd1 = [self.bin.bwa, "mem",
                "-t", str(num_threads),
                "-R", r'@RG\tID:ID_' + sample[0] + r'\tSM:' + sample[0],
                self.gref.bwa_db] + sample[1:3]

        cmd2 = [self.bin.samtools, "sort",
                "-o", out_bam_path]

        # Uncomment to skip if bam.bai already exists. For debuggin variant calling.
        #if os.path.isfile(out_bam_path + ".bai"):
        #    return(out_bam_path)

        proc1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        proc2 = subprocess.Popen(cmd2, stdin = proc1.stdout)
        proc2.communicate()
        return(out_bam_path)

    def rmdup_and_index(self, bam_path):
        # samptools rmdup and index for a given bam
        # Result is stored in same dir
        # Original bam file will be relaced

        bam_dir = os.path.dirname(bam_path)
        bam_name = os.path.basename(bam_path)

        tmp_bam_path = bam_dir + "/_" + bam_name + ".bam"

        cmd1 = [self.bin.samtools, "rmdup",
               "-s",
               bam_path,
               tmp_bam_path
               ]
        cmd2 = ["mv",
               tmp_bam_path,
               bam_path
               ]
        cmd3 = [self.bin.samtools, "index",
               bam_path
               ]

        subprocess.call(cmd1)
        subprocess.call(cmd2)
        subprocess.call(cmd3)

    def variant_analysis(self, out_dir, in_bam, sample_name):
        # Variant calling using gatk and annotion usig bcftools csq
        # Returns path of out_table

        out_gvcf = out_dir + "/" + sample_name + ".g.vcf"
        out_vcf = out_dir + "/" + sample_name + ".vcf"
        out_csqvcf = out_dir + "/" + sample_name + "_csq.vcf"
        out_table = out_dir + "/" + sample_name + "_out_table"

        cmd1 = [self.bin.gatk, "HaplotypeCaller",
               "-L", self.gref.bed,
               "-R", self.gref.ref_fa,
               "-I", in_bam,
               "-O", out_gvcf,
               "-ERC", "GVCF",
               "--max-mnp-distance", "5",
               "--native-pair-hmm-threads", "1"]

        cmd2 = [self.bin.gatk, "GenotypeGVCFs",
               "-R", self.gref.ref_fa,
               "-V", out_gvcf,
               "-O", out_vcf]

        cmd3 = [self.bin.bcftools, "csq",
               "-g", self.gref.gff3,
               "-f", self.gref.ref_fa,
               "-p", "a"
               "-l",
               "-o", out_csqvcf,
               out_vcf]

        proc1 = subprocess.call(cmd1)
        proc2 = subprocess.call(cmd2)
        proc3 = subprocess.call(cmd3)

        cmd4 = [self.bin.bcftools, "query",
               "-f", r"[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TBCSQ\t%TGT\t%AD\n]",
               "-o", out_table,
               out_csqvcf]

        subprocess.call(cmd4)
        return(out_table)

        #with open(final_out_table, 'w') as f:
        #    subprocess.Popen(cmd4, stdout = f)

import sys
import os

class GenomeRef:
    # This is a class to hold path information of each species
    def __init__(self, ref_dir, species):
        root_dir = os.path.join(ref_dir, species)
        self.bwa_db = os.path.join(root_dir, "bwadb", "ref")
        self.gff3 =  os.path.join(root_dir, "ref.gff3")
        self.ref_fa =  os.path.join(root_dir, "ref.fa")
        self.bed =  os.path.join(root_dir, "ref.bed")
        self.mdom_fa = os.path.join(root_dir, "ref.mdom.fa")
        self.mdom_coord_ck = os.path.join(root_dir, "ref.mdomcoord")
        self.mdom_coord_dl = os.path.join(root_dir, "ref.mdomcoord")

        self.check_existence()

    def check_existence(self):
        for file in [self.bwa_db + ".pac",
                     self.gff3,
                     self.ref_fa,
                     self.bed,
                     self.mdom_fa,
                     self.mdom_coord_ck]:
            if  os.path.exists(file):
                print("Using the reference file " + file + ".", file = sys.stderr)
            else:
                print(r"ERROR!: " + file + " does not exists /(*o*)\ !", file = sys.stderr)
                sys.exit(1)

class Bin:
    def __init__(self, script_dir):
        self.bwa = os.path.join(script_dir, "bin/bwa/bwa")
        self.samtools = os.path.join(script_dir, "bin/samtools-1.9/samtools")
        self.bcftools = os.path.join(script_dir, "bin/bcftools-1.9/bcftools")
        #self.freebayes = os.path.join(script_dir, "bin/freebayes/bin/freebayes")
        self.gatk = os.path.join(script_dir, "bin/gatk-4.0.11.0/gatk")

        self.check_existence()

    def check_existence(self):
        for key in self.__dict__.keys():
            file = self.__dict__[key]
            if os.path.isfile(file):
                print("Using " + key + " -> " + file, file = sys.stderr)
            else:
                print(r"ERROR!: " + file + " does not exists /(*o*)\ !", file = sys.stderr)
                sys.exit(1)

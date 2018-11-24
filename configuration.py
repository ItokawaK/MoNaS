import sys
import os
import subprocess
import shutil



class GenomeRef:
    '''
    This is a class to hold path information of references for each species.
    An instance is constructed with root directly of the MoNaS references and
    a code of species (name of the subdirectly) are given.
    The constructer automatically checks existence of bwa index and create them
    if absent.
    '''

    for program in ['bwa', 'samtools', 'bcftools', 'gatk']:
        if shutil.which(program) == None:
            print("Error: " + program + " is not in $PATH /(*o*)\ !")
            sys.exit(1)

    def __init__(self, ref_dir, species):
        root_dir = os.path.join(ref_dir, species)
        self.bwa_db = os.path.join(root_dir, "bwadb", "ref")
        self.gff3 =  os.path.join(root_dir, "ref.gff3")
        self.ref_fa =  os.path.join(root_dir, "ref.fa")
        self.bed =  os.path.join(root_dir, "ref.bed")
        self.mdom_fa = os.path.join(root_dir, "ref.mdom.fa")
        self.mdom_coord_ck = os.path.join(root_dir, "ref.mdomcoord")
        self.mdom_coord_dl = os.path.join(root_dir, "ref.mdomcoord")

        if not os.path.isfile(self.bwa_db  + ".pac"):
            print("Creating bwadb for " + self.ref_fa + "\n This may take for a while.", file = sys.stderr)
            self.creat_bwadb(root_dir)

        if not os.path.isfile(self.ref_fa + '.fai'):
            self.faidx()

        if not os.path.isfile(self.ref_fa + '.dict'):
            self.fa_dict()

        self.check_existence()


    def creat_bwadb(self, dir):
        bwadb = os.path.join(dir, 'bwadb')
        if not os.path.exists(bwadb):
            os.mkdir(bwadb)
        subprocess.call(['bwa', 'index',
                         '-p', bwadb + '/ref',
                        self.ref_fa])
    def faidx(self):
        subprocess.call(['samtools', 'faidx', self.ref_fa])

    def fa_dict(self):
        subprocess.call(['gatk', 'CreateSequenceDictionary',
                        '-R', self.ref_fa])

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

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

import sys
import os
import subprocess
import shutil
import json

from . import bed2gff3


class GenomeRef:
    '''
    This is a class to hold path information of references for each species.
    An instance is constructed with root directly of the MoNaS references and
    a code of species (name of the subdirectly) are given.
    The constructer automatically checks existence of bwa index and create them
    if absent.
    '''

    def __init__(self, ref_dir, species):

        self.root_dir = os.path.join(ref_dir, species)
        if not os.path.isdir(self.root_dir):
            print(self.root_dir + " was not found!", file = sys.stderr)
            sys.exit(1)

        self.bwa_db = os.path.join(self.root_dir, "bwadb", "ref")
        self.hisat_db = os.path.join(self.root_dir, "hisatdb", "ref")
        self.gff3 =  os.path.join(self.root_dir, "ref.gff3")
        self.ref_fa =  os.path.join(self.root_dir, "ref.fa")
        self.bed =  os.path.join(self.root_dir, "ref.bed")
        self.mdom_fa = os.path.join(self.root_dir, "ref.mdom.fa")

    #########################################################################
    def check_genomedb(self, mode, variant_caller, num_cpu = 1):
        if not os.path.isfile(self.ref_fa + '.fai'):
            subprocess.call(['samtools', 'faidx', self.ref_fa])

        #Creating a gff3 file if absent
        if not os.path.isfile(self.gff3):
            if os.path.isfile(self.bed):
                out_line = bed2gff3.create_bed_from_gff3(self.bed)
                with open(self.gff3, 'w') as f:
                    for l in out_line:
                        f.write(l + "\n")

        #Creating a protein alignment file if absent
        if not os.path.isfile(self.mdom_fa):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            cmd = [script_dir + "/make_AA_alignment.py",
                   "-o", self.mdom_fa,
                   self.ref_fa, self.bed]
            subprocess.call(cmd)

        if variant_caller == "gatk" and not os.path.isfile(os.path.join(self.root_dir, "ref.dict")):
            subprocess.call(['gatk', 'CreateSequenceDictionary',
                            '-R', self.ref_fa])

        if mode == 'ngs_dna':
            bwadb = os.path.join(self.root_dir, 'bwadb')
            if not os.path.exists(bwadb):
                print("Creating bwadb for " + self.ref_fa + "\n This may take for a while."\
                      " Please wait patiently.",
                      file = sys.stderr)
                os.mkdir(bwadb)
                subprocess.call(['bwa', 'index',
                                '-p', bwadb + '/ref',
                                self.ref_fa])
        elif mode == 'ngs_rna':
            hisatdb = os.path.join(self.root_dir, 'hisatdb')
            if not os.path.exists(hisatdb):
                print("Creating hisat2db for " + self.ref_fa + "\n This may take for a while."\
                      " Please wait patiently.",
                      file = sys.stderr)
                os.mkdir(hisatdb)
                subprocess.call(['hisat2-build',
                                '-p', str(num_cpu),
                                self.ref_fa,
                                hisatdb + '/ref'])
        else:
            print(mode + " is currently not supproted.")
            sys.exit(1)

    def check_program_path(self, mode, variant_caller):

        with open(os.path.dirname(os.path.abspath(__file__)) + '/bin_path.json') as f:
            program_path = json.load(f)

        required_programs = ['samtools', 'bcftools']

        if mode == 'ngs_dna':
            required_programs.append('bwa')
        if mode == 'ngs_rna':
            required_programs.append('hisat2')

        required_programs.append(variant_caller)

        for program in required_programs:
            if (program in program_path) and program_path[program]:
                abs_path = os.path.expanduser(program_path[program])
                if os.path.isfile(os.path.join(abs_path, program)):
                    os.environ['PATH'] = abs_path + os.pathsep + os.environ['PATH']
                else:
                    print(os.path.join(abs_path, program) + " was not found /(*o*)\ !",
                          file = sys.stderr)
                    sys.exit(1)
            else:
                if shutil.which(program) == None:
                    print("Error: " + program + " is not in $PATH /(*o*)\ !",
                          file = sys.stderr)
                    sys.exit(1)

    def check_existence_for_genotype(self):
        for file in [self.gff3,
                     self.ref_fa,
                     self.bed,
                     self.mdom_fa]:
            if  os.path.exists(file):
                print("Using the reference file " + file + ".", file = sys.stderr)
            else:
                print(r"ERROR!: " + file + " does not exists /(*o*)\ !", file = sys.stderr)
                sys.exit(1)

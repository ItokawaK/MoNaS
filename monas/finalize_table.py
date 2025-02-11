#!/usr/bin/env python

import subprocess as sp
import re
import os
import sys
import json

'''
Copyright (c) 2024, Kentaro Itokawa <itokawa@niid.go.jp>

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

class TBCSQ():
    '''
    This object recieve translated BCSQ string from BCFtools query for one allele, and MDom_convert instance 
    to convert aa position of arbitrary insect to M. domestica
    functions:
      @property aa_change_str(): return position unconverted aa change notion like F1550C
      @property aa_change_str_mdom(): return position converted aa change notion like F1534C
    '''
    
    def __init__(self, tbcsq_line, mdom_convert):
        self.raw = tbcsq_line

        self.type = 'wild' # e.g. missense, stop_gained, frameshift, etc
        self.aa_pos = None # aa position in VCF
        self.aa_ref = None # reference aa residue
        self.aa_alt = None # mutated aa residue
        
        self.aa_pos_mdom = None # converted aa position by MDom_convert instance 
        self.kdr_evidence = 'NA' # from MDom_convert instance 

        self.cds_change = False # True if the mutaion affect CDS nucleotide 
        
        if tbcsq_line == '.': # allele is wildtype nucleotide
            return
            
        self.csq_els = tbcsq_line.split(',') 
        for csq_el in self.csq_els:
            els = csq_el.split('|')
            type = els[0] 
            if type in ['missense', 'stop_gained', 'frameshift']: # aa changing mutation
                aa_change_str = els[5] # e.g. 1234M>1234D
                matchOB = re.match(r"(\d+)([^\d]+)>(\d+)([^\d]+)", aa_change_str)
                if matchOB:
                    self.aa_pos = int(matchOB.group(1))
                    self.aa_ref = matchOB.group(2)
                    self.aa_alt = matchOB.group(4)

                mdom_pos, kdr_evidence = mdom_convert.mos2fly(self.aa_ref, self.aa_alt, str(self.aa_pos))
                self.aa_pos_mdom = mdom_pos
                self.kdr_evidence = kdr_evidence
                self.cds_change = True
                self.type = type
                # self.aa_change = aa_change
                break
                
            if type == 'synonymous':
                aa_change_str = els[5] # e.g. 1234M
                matchOB = re.match(r"(\d+)([^\d]+)", aa_change_str) 
                self.aa_pos = int(matchOB.group(1))
                self.aa_ref = self.aa_alt = matchOB.group(2)
                mdom_pos, kdr_evidence = mdom_convert.mos2fly(self.aa_ref, self.aa_alt, str(self.aa_pos))
                self.aa_pos_mdom = mdom_pos
                self.cds_change = True
                self.type = type
                break

    @property
    def aa_change_str(self):
        if  self.cds_change:
            if self.type != 'synonymous':
                return f'{self.aa_ref}{self.aa_pos}{self.aa_alt}'
            else:
                return f'{self.aa_pos}{self.aa_alt}'
        else:
            return 'wild'

    @property
    def aa_change_str_mdom(self):
        if  self.cds_change:
            if self.type != 'synonymous':
                return f'{self.aa_ref}{self.aa_pos_mdom}{self.aa_alt}'
            else:
                return f'Synonymous'
        else:
            return 'wild'
        
class BCF_QUERY():
    '''
    This object recieve a string generated from bcftools query function, MDom_convert 
    instacen and Bed instance
    Upon initialization, it calculate aa position conversion, exon position
    functions:
      @property tbl_line(): Output tab delimitted line for final table
    '''
    
    _queries_ = '%SAMPLE %CHROM %POS %REF %ALT %QUAL %TGT %TBCSQ %AD'.split(' ')
    query_str = "[" + '\t'.join(_queries_) + "\n]"

    tbl_headr = ["#ID",
                  "CHROM",
                  "POS",
                  "REF_ALLELE",
                  "ALT_ALLELE",
                  "QUAL",
                  "GENOTYPE",
                  "AA_CHANGE",
                  "AA_CHANGE_MDOM",
                  "KDR_EVIDENCE",
                  "AD",
                  "EXON"]
    
    def __init__(self, in_line, mdom_convert, bed):
        els = in_line.split('\t')
        self.sample = els[0] 
        self.chrom = els[1]
        self.pos = int(els[2])
        self.nuc_ref = els[3]
        self.nuc_alt = els[4]
        self.qual = els[5]
        self.gt = els[6]
        self.tbcsq_al1 = TBCSQ(els[7], mdom_convert)
        self.tbcsq_al2 = TBCSQ(els[8], mdom_convert)
        self.ad = els[9]
        
        self.aa_gt = self.tbcsq_al1.aa_change_str + '/' + self.tbcsq_al2.aa_change_str # aa genotype, position unconverted
        self.type_gt =  self.tbcsq_al1.type + '/' + self.tbcsq_al2.type # mutation type genotype

        self.aa_gt_mdom =  self.tbcsq_al1.aa_change_str_mdom + '/' + self.tbcsq_al2.aa_change_str_mdom # aa genotype, position converted
        self.kdr_evidence_gt =  self.tbcsq_al1.kdr_evidence + '/' + self.tbcsq_al2.kdr_evidence # kdr evidence genotype
        
        self.exon = bed.which_exon(self.pos) # corresponding exon
        
    @property
    def tbl_line(self):
        out_fields = [self.sample,          # 1
                      self.chrom,           # 2 
                      str(self.pos),        # 3
                      self.nuc_ref,         # 4
                      self.nuc_alt,         # 5
                      self.qual,            # 6
                      self.gt,              # 7
                      self.aa_gt,           # 8
                      self.aa_gt_mdom,      # 9
                      self.kdr_evidence_gt, # 10
                      self.ad,              # 11
                      self.exon]            # 12
        
        return '\t'.join(out_fields)
    
class Bed:    
    '''
    This object represent information in BED file for CDS structure
    functions:
      which_exon(pos) # ask which exon the pos nucleotide is included
    '''
    class Exon:
        def __init__(self, start , end, name, strand):
            self.name = name
            self.no = int(re.search(pattern= "[0-9]+", string= name).group())
            self.start = int(start)
            self.end = int(end)
            self.strand = strand
            if name[-1] in ('c', 'k', 'd', 'l'):
               self.mut_excl = name[-1]
            else:
               self.mut_excl = None

    def __init__(self, bed_file):

        with open(bed_file) as f:
            bed_lines = [l.rstrip() for l in f.readlines() if not l.startswith("#")]

        self.entries = []
        self.headers = []
        for line in bed_lines:
            last_strand = None
            if line.startswith('browser') or  line.startswith('track'):
                self.headers.append(line.split("\t"))
            else:
                tmp = line.split("\t")
                self.chrom = tmp[0]
                if len(tmp) > 5:
                    self.entries.append(self.Exon(tmp[1], tmp[2], tmp[3], tmp[5]))
                    fix_strand = False
                    if last_strand and last_strand != tmp[5]:
                        sys.exit("Error!: There are two different strands in " + bed_file)
                    last_strand = tmp[5]
                else:
                    self.entries.append(self.Exon(tmp[1], tmp[2], tmp[3], "+"))
                    fix_strand = True

        if fix_strand:
            if (self.entries[0].no < self.entries[1].no) == \
                                (self.entries[0].start < self.entries[1].end):
                self.strand = "+"
            else:
                self.strand = "-"
            for exon in self.entries:
                exon.strand = self.strand
        else:
            self.strand = self.entries[0].strand

    def which_exon(self, pos):
            for Exon in self.entries:
                if pos > Exon.start and pos <= Exon.end:
                    return(Exon.name)
            return('intron')
    
class MDom_convert:
    '''
    This object recieve multi alingment fasta file for vgsc protein of M.domesticat and 
    species of your interest. Path of kdr_evidence.json file.
    functions:
      mos2fly(aa_ref, aa_alt, aa_pos) # Return tupple ("converted aa position"[int], "kdr evidence"[str])
    
    '''

    
    def __init__(self, fasta, kdr_list_path): #kdr_list_path is a json file
        self.hash = {} # key: unconverted position (str), value: converted position (str) "-" if no corresponding aa
        self._fasta_2_hash_(fasta)

        with open(kdr_list_path) as f:
            self.kdr_dict = json.load(f)

    def _fasta_2_hash_(self, fasta):
        """
        Creating [mosquito AA pos] -> [M. domestica AA pos] hash function from
        multiple aligment fasta. The first entry in the fasta will be regarded
        as M. domestica AA seq. Currently, the third and follwing sequences will
        be ignored because it assumes no indels between transcript variants 
        """

        with open(fasta) as f:
            AA_seqs = []

            if f.readline().startswith(">"):
                seq_buf = []
            else:
                print(fasta + ' is not in valid format!', file = sys.stderr)
                sys.exit(1)

            for line in f.readlines():
                if line.startswith(">"):
                    AA_seqs.append(seq_buf)
                    seq_buf = []
                else:
                    seq_buf += list(line.rstrip())

            AA_seqs.append(seq_buf)



        if not len(AA_seqs[0]) == len(AA_seqs[1]):
            raise ValueError("error!")

        idx = [0, 0]

        for mdom, mos in zip(AA_seqs[0], AA_seqs[1]):
            if mos != "-":
                idx[0] += 1
                mos_no = str(idx[0])
            else:
                mos_no = "-"
            if mdom != "-":
                idx[1] += 1
                mdom_no = str(idx[1])
            else:
                mdom_no = "-"

            if mos != "-":
                self.hash[mos_no] = mdom_no

            # if debug:
            #     print("\t".join([str(idx[0]), mos, mdom, str(idx[1])]), file = sys.stderr)

    def mos2fly(self, aa_ref, aa_alt, aa_pos):
        kdr_evidence = "NA"

        # if AA_str == "wild":
        #     return ("wild", kdr_evidence)

        if aa_ref == aa_alt: # if it is synonymous 
            kdr_evidence = 'NA'
            mdom_aa_pos = self.hash[aa_pos]
            return (mdom_aa_pos, kdr_evidence)

        # otherwise
        mdom_aa_pos = self.hash[aa_pos]
        kdr_evidence = "Unknown"
        if mdom_aa_pos in self.kdr_dict: 
            kdr_evidence = "Presumable" # meaning the residue is pyR but the effect of alt aa is not known 
            if aa_alt in self.kdr_dict[mdom_aa_pos]:
                kdr_evidence = self.kdr_dict[mdom_aa_pos][aa_alt]

        return (mdom_aa_pos, kdr_evidence)
    

def crate_table(vcf_path, bed_file, mdom_aln_fasta, out_table_file):

    # Run bcftools query    
    proc = sp.run(['bcftools', 'query', '-f', BCF_QUERY.query_str, vcf_path], stdout=sp.PIPE)

    # kdr_list.json path, which should be under the same dir to this script
    kdr_list = os.path.join(os.path.dirname(os.path.abspath(__file__)) , "kdr_list.json")

    md_conv = MDom_convert(mdom_aln_fasta, kdr_list)
    bed = Bed(bed_file)

    # entries in table header
    header = BCF_QUERY.tbl_headr

    with open(out_table_file, 'w') as f:
        print('\t'.join(header), file=f)
        for l in proc.stdout.decode().rstrip().split('\n'):
            bcf_tools_q = BCF_QUERY(l, md_conv, bed)
            print(bcf_tools_q.tbl_line, file=f)

def main(args):
    crate_table(args.out_csq, args.ref_bed, args.mdom_fa, "/dev/stdout")
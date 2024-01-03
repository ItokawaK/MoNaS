#!/usr/bin/env python

import os
import sys
import re
import json
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

debug = False

class VCF_line:
    """
    This object represents a single line of VCF file with TBSQ field.
    For initialization, follwing argments are required.
       line: a string of VCF lines
       bed: BED object created from ref.bed file
       mdom_conv: a MDom_comvert object created from ref.mdom.fa
       samples: a list of all sample names
    """
    def __init__(self, line, bed, mdom_conv, samples):

        def parse_GT(Fs, vcf_format, samples):
            """
            Fuction to create dictionary of dictionary to represent
            per sample GT information.
            Arguments:
              Fs: list of genotype fields
              vcf_format: list of FORMAT elements (GT, AD, DP, etc.)
              samples: list of sample name in same order to VCF file
            Output:
              dictionary:
                  key: sample name, value:dictionary2
                  dictionary2:
                     key: FORMAT element, value: value in GT field.
            """
            out_dict = {}
            for i in range(len(samples)):
                sample_dict = {}
                sample_F = Fs[i].split(":")
                for ii in range(len(vcf_format)):
                    sample_dict[vcf_format[ii]] = sample_F[ii]
                out_dict[samples[i]] = sample_dict
            return out_dict
        self.line = line
        Fs = self.line.split("\t")

        self.entries = dict()
        self.entries['CHROM'] = Fs[0]
        self.entries['POS'] = int(Fs[1])
        self.entries['REF_ALLELE'] = Fs[3]
        self.entries['ALT_ALLELE'] = Fs[4]
        self.entries['QUAL'] = float(Fs[5])
        self.entries['EXON'] = bed.which_exon(self.entries['POS'])

         # list of all alleles on this position
        self.alleles = [self.entries['REF_ALLELE']] + self.entries['ALT_ALLELE'].split(",")
        # dictionary representing INFO field values
        info = {ll[0]:ll[1] for ll in [l.split("=") for l in Fs[7].split(";")]}
        if "BCSQ" in info:
            self.tbcsq = BCSQ(info["BCSQ"])
        else:
            self.tbcsq = None
        #self.exon = bed.which_exon(self.entries['POS'])
        self.mdom_conv = mdom_conv
        format = Fs[8].split(":")
        self.sample_data = parse_GT(Fs[9:], format, samples)
        for sample in samples:
            self.set_sample_data(sample)

    def set_sample_data(self, sample_name):
        """
        This method is used to parse BCSQ field and GT, determining AA genotype
        of each samples. Set per sample data (GENOTYPE, AA_CHANGE, etc.)
        """

        if self.sample_data[sample_name]["GT"] == ".":
            return None
        else:
            GT_indx = [int(i) for i in self.sample_data[sample_name]["GT"].split("/")]

        alleles = ["", ""] # alleles in individual genotype
        AA_change = ["", ""] # corresponding AA changes (e.g. 123F>123L) to each allele
        AA_change_simple = ["", ""] # corresponding AA changes (e.g. F123L) to each allele
        AA_change_mdom  = ["", ""] # corresponding Mdom AA changes (e.g. F120L) to each allele
        kdr_type =  ["", ""] # corresponding kdr type (Unknown, Strong, Supportive) to each AA change
        for i in [0, 1]:
            alleles[i] = self.alleles[GT_indx[i]]

            if alleles[i] == self.entries['REF_ALLELE']:
                AA_change[i] = "wild"
                continue

            if len(self.tbcsq.dict_ck) == 0 and len(self.tbcsq.dict_dl) == 0:
                AA_change[i] = "synonymous"
            elif len(self.tbcsq.dict_ck) > 0:
                try:
                    AA_change[i] = self.tbcsq.dict_ck[alleles[i]]
                except:
                    print("error!:" + self.line)

            else:
                try:
                    AA_change[i] = self.tbcsq.dict_dl[alleles[i]]
                except:
                    print("error!:" + self.line)


        for i in [0, 1]:
            matchOB = re.match(r"(\d+)([^\d]+)>(\d+)([^\d]+)", AA_change[i])
            if matchOB:
                coord = matchOB.group(1)
                old_AA = matchOB.group(2)
                new_AA = matchOB.group(4)
                AA_change_simple[i] = old_AA + coord + new_AA
            else:
                AA_change_simple[i] = AA_change[i]

            AA_change_mdom[i], kdr_type[i] = self.mdom_conv.mos2fly(AA_change[i])

        self.sample_data[sample_name]["GENOTYPE"] = alleles[0] + "/" + alleles[1]
        self.sample_data[sample_name]["AA_CHANGE"] = "/".join(AA_change_simple)
        self.sample_data[sample_name]["AA_CHANGE_MDOM"] = "/".join(AA_change_mdom)
        self.sample_data[sample_name]["KDR_EVIDENCE"] = "/".join(kdr_type)

    def get_info(self, sample, info_to_get):
        """
        Get information in info_to_get for sample
        info_to_get is a list of string formated like
         ["CHROM", "POS", "REF_ALLELE", ...]
        """
        out_list = []
        for info in info_to_get:
            if info in self.entries:
                out_list.append(self.entries[info])
            elif info in self.sample_data[sample]:
                out_list.append(self.sample_data[sample][info])
            else:
                out_list.append("NO_DATA")
        return out_list

class BCSQ:
    """
    A class for INFO/BCSQ field of a single VCF line
    """

    def __init__(self, info_str):

        self.dict_ck = {} # key: ALT_allele, value: AA subsutuitoin eg. nnnX>nnnX
        self.dict_dl = {} # key: ALT_allele, value: AA subsutuitoin eg. nnnX>nnnX
        #self.allele = 'Wild'

        if info_str != '.':
            split_info = info_str.split(",")
            split_info_lists = [tmp.split("|") for tmp in split_info]

            for info_list in split_info_lists:
                do_skip = False
                for pat in ["splice_region", "splice_acceptor", "splice_donor", "start_lost"]:
                    if pat in info_list[0]:
                        do_skip = True
                if do_skip:
                    continue
                if(info_list[0].startswith("@")):
                    continue
                splice_variant = info_list[2][-2:] # "ck" or "dl"
                try:
                    alt_allele = info_list[6].split(">")[1]
                except:
                    print("Unexpexted format of vcf: " + info_str, file = sys.stderr)
                    print(info_list, file = sys.stderr)

                if splice_variant == "ck":
                    self.dict_ck[alt_allele] = info_list[5]

                if splice_variant == "dl":
                    self.dict_dl[alt_allele] = info_list[5]


class MDom_comvert:
    def __init__(self, fasta, kdr_list_path): #kdr_list_path is a json file
        self.hash = {}
        self.fasta_2_hash(fasta)

        with open(kdr_list_path) as f:
            self.kdr_dict = json.load(f)

    def fasta_2_hash(self, fasta):
        """
        Creating [mosquito AA pos] -> [M. domestica AA pos] hash function from
        multiple aligment fasta. The first entry in the fasta will be regarded
        as M. domestica AA seq. Currently, the third and follwing sequences will
        be ignored because ck dl transcript variants will have exactly same
        aligment to M. domestica VGSC.
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

            if debug:
                print("\t".join([str(idx[0]), mos, mdom, str(idx[1])]), file = sys.stderr)

    def mos2fly(self, AA_str):
        kdr_evidence = "NA"

        if AA_str == "wild":
            return ("wild", kdr_evidence)

        match_change = re.match(r"(\d+)([^\d]+)>(\d+)([^\d]+)", AA_str)
        if match_change:
            mos_coord = match_change.group(1)
            old_AA = match_change.group(2)
            new_AA = match_change.group(4)
            fly_coord = self.hash[mos_coord]

            if not (fly_coord + new_AA) in self.kdr_dict:
                kdr_evidence = "Unknown"
            elif self.kdr_dict[(fly_coord + new_AA)] == 1:
                kdr_evidence = "Strong"
            elif self.kdr_dict[(fly_coord + new_AA)] == 2:
                kdr_evidence = "Supportive"
            return (old_AA + fly_coord + new_AA, kdr_evidence)
        else:
            return ("synonymous", kdr_evidence)


class Bed:
    class Exon:
        def __init__(self, start , end, name, strand):
            self.name = name
            self.no = int(re.search(pattern= "\d+", string= name).group())
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

def create_table(csqvcfs, info_to_get, bed_file, mdom_fasta, out_table_file, header = 'auto'):
    """
    Out put table in single out_table_file with header containig info_to_get

    1. prints header string
    2. Reads VCF files in list csqvcfs
    3. Retrieves sample name from line starts with #CHROM
    4. Converts VCF lines into VCF_line object
    5. Retrieves information described in info_to_get string from VCF_line object
    """

    def vcf_to_table(in_fh, bed, md_conv, info_to_get):
        out_lines = []
        for l in in_fh.readlines():
            if l.startswith("#CHROM"):
                samples = l.rstrip().split("\t")[9:]
            if l.startswith("#"):
                continue

            vcf_l = VCF_line(l.rstrip(), bed, md_conv, samples)
            if vcf_l.tbcsq == None:
                continue
            for sample in samples:
                if vcf_l.sample_data[sample]["GT"] in ["0/0", "."]:
                    continue

                info = vcf_l.get_info(sample, info_to_get)
                out_str = [sample] + [str(info_) for info_ in info]
                out_lines.append("\t".join(out_str))
        return out_lines

    kdr_list = os.path.dirname(os.path.abspath(__file__)) + "/kdr_list.json"
    md_conv = MDom_comvert(mdom_fasta, kdr_list)
    bed = Bed(bed_file)

    if header == 'auto':
        header = "\t".join(["#ID"] + info_to_get)

    print_header = True

    with open(out_table_file, "w") as out_f:
        for csqvcf in csqvcfs:
            with open(csqvcf) as f:
                if print_header:
                    print(header, file = out_f)
                print_header = False
                write_us = vcf_to_table(f, bed, md_conv, info_to_get)
                for l in write_us:
                    print(l, file = out_f)

def main():
    import argparse

    parser = argparse.ArgumentParser(description = 'Convert csqvcf to table')

    parser.add_argument("out_csq",
                        help = "out_csq file path")

    parser.add_argument("ref_bed",
                        help = "ref.bed file path")

    parser.add_argument("mdom_fa",
                        help = "ref.mdom.fa file path")

    args = parser.parse_args()

    create_table(csqvcfs = [args.out_csq],
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
                 bed_file = args.ref_bed,
                 mdom_fasta = args.mdom_fa,
                 out_table_file = "/dev/stdout"
                 )
    
if __name__ == '__main__':
    main()
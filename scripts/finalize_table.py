#!/usr/bin/env python

import os
import sys
import re
import json

debug = False

class VCF_line:
    def __init__(self, line, bed, mdom_conv, samples):

        def vcf_samples(Fs, vcf_format, samples):
            out_dict = {}
            for i in range(len(samples)):
                sample_dict = {}
                sample_F = Fs[i].split(":")
                for ii in range(len(vcf_format)):
                    sample_dict[vcf_format[ii]] = sample_F[ii]
                out_dict[samples[i]] = sample_dict
            return out_dict

        Fs = line.split("\t")

        self.entries = dict()
        self.entries['CHROM'] = Fs[0]
        self.entries['POS'] = int(Fs[1])
        self.entries['REF_ALLELE'] = Fs[3]
        self.entries['ALT_ALLELE'] = Fs[4]
        self.entries['QUAL'] = float(Fs[5])
        self.entries['EXON'] = bed.which_exon(self.entries['POS'])


        self.line = line
        self.alleles = [Fs[3]] + Fs[4].split(",")
        info = {ll[0]:ll[1] for ll in [l.split("=") for l in Fs[7].split(";")]}
        if "BCSQ" in info:
            self.tbcsq = BCSQ(info["BCSQ"])
        else:
            self.tbcsq = None
        self.exon = bed.which_exon(self.entries['POS'])
        self.mdom_conv = mdom_conv
        format = Fs[8].split(":")
        self.sample_data = vcf_samples(Fs[9:], format, samples)
        for sample in samples:
            self.set_sample_data(sample)

    def set_sample_data(self, sample_name):

        if self.sample_data[sample_name]["GT"] == ".":
            return None

        GT_indx = [int(i) for i in self.sample_data[sample_name]["GT"].split("/")]
        al = [ "", ""]
        AA_change = ["", ""]
        for i in [0, 1]:
            al[i] = self.alleles[GT_indx[i]]

            if al[i] == self.entries['REF_ALLELE']:
                AA_change[i] = "wild"
                continue

            if len(self.tbcsq.dict_ck) == 0 and len(self.tbcsq.dict_ck) == 0:
                AA_change[i] = "synonymous"
            elif len(self.tbcsq.dict_ck) > 0:
                try:
                    AA_change[i] = self.tbcsq.dict_ck[al[i]]
                except:
                    print("error!:" + self.line)

            else:
                try:
                    AA_change[i] = self.tbcsq.dict_dl[al[i]]
                except:
                    print("error!:" + self.line)

        AA_change_mdom = [self.mdom_conv.mos2fly(AA_change[0]),
                          self.mdom_conv.mos2fly(AA_change[1])]

        self.sample_data[sample_name]["GT"] = al[0] + "/" + al[1]
        self.sample_data[sample_name]["AA_CHANGE"] = \
                                    AA_change[0] + "/" + AA_change[1]
        self.sample_data[sample_name]["AA_CHANGE_MDOM"] = \
                                  AA_change_mdom[0] + "/" + AA_change_mdom[1]

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
                for pat in ["splice_region", "splice_acceptor", "splice_donor"]:
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
        Mdom_AA_aligned = []
        Mos_AA_aligned = []

        with open(fasta) as f:
            if f.readline().startswith(">"):
                is_Mdom = True
            else:
                print(fasta + ' is not in valid format!', file = sys.stderr)
                sys.exit(1)
            for line in f.readlines():
                if line.startswith(">"):
                    is_Mdom = False
                    continue
                if is_Mdom:
                    Mdom_AA_aligned += list(line.rstrip())
                else:
                    Mos_AA_aligned += list(line.rstrip())

        if not len(Mdom_AA_aligned) == len(Mos_AA_aligned):
            raise ValueError("error!")

        idx = [0, 0]

        for mos, mdom in zip(Mos_AA_aligned, Mdom_AA_aligned):
            if mos != "-":
                idx[0] += 1
                a = str(idx[0])
            else:
                a = "-"
            if mdom != "-":
                idx[1] += 1
                b = str(idx[1])
            else:
                b = "-"

            if a != "-":
                self.hash[a] = b

            if debug:
                print("\t".join([str(idx[0]), mos, mdom, str(idx[1])]), file = sys.stderr)

    def mos2fly(self, AA_str):
        if AA_str == "wild":
            return("wild")

        matchOB = re.match(r"(\d+)([^\d]+)>(\d+)([^\d]+)", AA_str)
        if matchOB:
            mos_coord = matchOB.group(1)
            old_AA = matchOB.group(2)
            new_AA = matchOB.group(4)
            fly_coord = self.hash[mos_coord]
            kdr_symbol = ""
            if (fly_coord + new_AA) in self.kdr_dict:
                if self.kdr_dict[(fly_coord + new_AA)] == 1:
                    kdr_symbol = "!!"
                else:
                    kdr_symbol = "??"
            return(old_AA + fly_coord + new_AA + kdr_symbol)
        else:
            return("synonymous")


class Bed:
    class Exon:
        def __init__(self, name, start , end):
            self.name = name
            self.no = int(re.search(pattern= "\d+", string= name).group())
            self.start = int(start)
            self.end = int(end)
            if name[-1] in ('c', 'k', 'd', 'l'):
               self.mut_excl = name[-1]
            else:
               self.mut_excl = None

    def __init__(self, bed_file):

        with open(bed_file) as f:
            bed_lines = [l.rstrip() for l in f.readlines() if not l.startswith("#")]

        entries = []
        for line in bed_lines:
            self.chrom, start, end, name = line.split("\t")
            entries.append(self.Exon(name, start, end))
        self.entries = entries
        if (self.entries[1].no < self.entries[2].no) == \
                            (self.entries[1].start < self.entries[2].end):
            self.strand = "+"
        else:
            self.strand = "-"

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

    kdr_list = os.path.dirname(os.path.abspath(__file__)) + "/kdr_list.json"
    #print(kdr_list, file = sys.stderr)
    md_conv = MDom_comvert(mdom_fasta, kdr_list)
    bed = Bed(bed_file)

    if header == 'auto':
        header = "\t".join(["#ID"] + info_to_get)

    with open(out_table_file, "w") as out_f:
        print_header = True
        for csqvcf in csqvcfs:
            with open(csqvcf) as f:
                if print_header:
                    print(header, file = out_f)
                print_header = False
                for l in f.readlines():
                    if l.startswith("#CHROM"):
                        samples = l.rstrip().split("\t")[9:]
                    if not l.startswith("#"):
                        vcf_l = VCF_line(l.rstrip(), bed, md_conv, samples)
                        if vcf_l.tbcsq == None:
                            continue
                        for sample in samples:
                            if not vcf_l.sample_data[sample]["GT"] in ["0/0", "."]:
                                out_str = [sample] + [str(info) for info in vcf_l.get_info(sample, info_to_get)]

                                print("\t".join(out_str), file = out_f)

if __name__ == '__main__':
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
                                "GT",
                                "AA_CHANGE",
                                "AA_CHANGE_MDOM",
                                "AD",
                                "EXON"],
                 bed_file = args.ref_bed,
                 mdom_fasta = args.mdom_fa,
                 out_table_file = "/dev/stdout"
                 )

#!/usr/bin/env python

import os
import sys
import re

debug = True

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
        self.line = line
        self.chrom = Fs[0]
        self.pos = int(Fs[1])
        self.ref = Fs[3]
        self.alt = Fs[4]
        self.alleles = [Fs[3]] + Fs[4].split(",")
        self.qual = float(Fs[5])
        self.info = {ll[0]:ll[1] for ll in [l.split("=") for l in Fs[7].split(";")]}
        if "BCSQ" in self.info:
            self.tbcsq = BCSQ(self.info["BCSQ"])
        else:
            self.tbcsq = None
        self.format = Fs[8].split(":")
        self.exon = bed.which_exon(self.pos)
        self.mdom_conv = mdom_conv
        self.sample_data = vcf_samples(Fs[9:], self.format, samples)


    def get_sample_data(self, sample_name):

        if self.sample_data[sample_name]["GT"] == ".":
            return None

        GT_indx = [int(i) for i in self.sample_data[sample_name]["GT"].split("/")]
        al = [ "", ""]
        AA_change = ["", ""]
        for i in [0, 1]:
            al[i] = self.alleles[GT_indx[i]]

            if al[i] == self.ref:
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


        return [sample_name,
                self.chrom,
                str(self.pos),
                self.ref,
                self.alt,
                al[0] + "/" + al[1],
                AA_change[0] + "/" + AA_change[1],
                self.mdom_conv.mos2fly(AA_change[0]) + "/" + self.mdom_conv.mos2fly(AA_change[1]),
                self.sample_data[sample_name]["AD"],
                self.exon
                 ]

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
    def __init__(self, fasta):
        self.hash = {}
        self.fasta_2_hash(fasta)

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
            return(old_AA + fly_coord + new_AA)
        else:
            return("synonymous")


class Bed:
    class Exon:
        def __init__(self, name, start , end):
            self.name = name
            self.start = int(start)
            self.end = int(end)

    def __init__(self, bed_file):

        with open(bed_file) as f:
            bed_lines = [l.rstrip() for l in f.readlines()]

        entries = []
        for line in bed_lines:
            chrom, start, end, name = line.split("\t")
            entries.append(self.Exon(name, start, end))
        self.entries = entries

    def which_exon(self, pos):
            for Exon in self.entries:
                if pos > Exon.start and pos <= Exon.end:
                    return(Exon.name)
            return('intron')

def create_table(csqvcf, bed_file, fasta, out_table_file):
    md_conv = MDom_comvert(fasta)
    bed = Bed(bed_file)

    with open(csqvcf) as f:
        with open(out_table_file, "w") as out_f:
            for l in f.readlines():
                if l.startswith("#CHROM"):
                    samples = l.rstrip().split("\t")[9:]
                if not l.startswith("#"):
                    vcf_l = VCF_line(l.rstrip(), bed, md_conv, samples)
                    if vcf_l.tbcsq == None:
                        continue
                    for s in samples:
                        if not vcf_l.sample_data[s]["GT"] in ["0/0", "."]:
                            print("\t".join(vcf_l.get_sample_data(s)), file = out_f)

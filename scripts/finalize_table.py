#!/usr/bin/env python

import os
import sys
import re

class VCF_info:

    def __init__(self, info_str):

        self.dict = {}
        self.allele = 'Wild'

        if info_str != '.':
            split_info = info_str.split(",")
            split_info_lists = [tmp.split("|") for tmp in split_info]

            for info_list in split_info_lists:
                if(info_list[0] in ["splice_region", "splice_acceptor", "splice_donor"]):
                    continue
                if(info_list[0].startswith("@")):
                    continue
                splice_variant = info_list[2][-2:]
                self.dict[splice_variant] = info_list
                self.allele = info_list


class VCF_line:

    def __init__(self, line, bed, mdom_coord):

        Fs = line.split("\t")
        #print(Fs)
        self.sample_id = Fs[0]
        self.chrom = Fs[1]
        self.pos = Fs[2]
        self.ref_allele = Fs[3]
        self.alt_alleles = Fs[4].split(",")
        self.qual = Fs[5]
        self.info = []
        self.info.append(VCF_info(Fs[6]))
        self.info.append(VCF_info(Fs[7]))
        self.gt = Fs[8].split("/")
        self.depth  = Fs[9]
        self.exon = bed.which_exon(self.pos)
        self.mdom_coord = mdom_coord

    def get_table(self):

        allele = []
        mdom_allele = []
        for i in [0, 1]:
            if self.info[i].allele == 'Wild':
                al = "Wild"
                md_al = "Wild"
            else:
                al = self.info[i].allele[5]
                md_al = self.mdom_coord.mos2fly(al)
            allele.append(al)
            mdom_allele.append(md_al)
        return("\t".join(
               [
                self.sample_id,
                self.chrom,
                self.pos,
                self.ref_allele,
                ",".join(self.alt_alleles),
                self.qual,
                "/".join(self.gt),
                allele[0] + "/" + allele[1],
                mdom_allele[0]  + "/" + mdom_allele[1],
                self.exon,
                self.depth
                ]
             )
        )

class Exon:
    def __init__(self, name, start , end):
        self.name = name
        self.start = start
        self.end = end

class Bed:
    def __init__(self, bed_file):

        with open(bed_file) as f:
            bed_lines = [l.rstrip() for l in f.readlines()]

        entries = []
        for line in bed_lines:
            chrom, start, end, name = line.split("\t")
            entries.append(Exon(name, start, end))
        self.entries = entries

    def which_exon(self, pos):
            for Exon in self.entries:
                if pos > Exon.start and pos <= Exon.end:
                    return(Exon.name)
            return('intron')


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
                    Mdom_AA_aligned += list(line)
                else:
                    Mos_AA_aligned += list(line)

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

    def mos2fly(self, AA_str):
        matchOB = re.match(r"(\d+)([^\d]+)>(\d+)([^\d]+)", AA_str)
        if matchOB:
            mos_coord = matchOB.group(1)
            old_AA = matchOB.group(2)
            new_AA = matchOB.group(4)
            fly_coord = self.hash[mos_coord]
            return(old_AA + fly_coord + new_AA)
        else:
            return("synonymous")

def create_table(table_file, bed_file, fasta, out_table_file):
    with open(table_file) as f:
        lines = f.readlines()

    bed = Bed(bed_file)
    mdom_coord = MDom_comvert(fasta)
    lines = [VCF_line(l.rstrip(), bed, mdom_coord) for l in lines]

    with open(out_table_file, 'w') as f:
        for l in lines:
            f.write(l.get_table() + "\n")

#!/usr/bin/env python

import sys
import os
from finalize_table import Bed
import argparse

def create_bed_from_gff3(bed_file):
    bed = Bed(bed_file)

    exons = bed.entries

    if bed.strand == "-":
        exons.sort(key = lambda ex: ex.start, reverse= True)

    def phase(n):
         return (3 - (n % 3)) % 3

    total_l = 0
    for exon in exons:
        if not exon.mut_excl in ('d', 'l'):
            exon.phase = phase(total_l)
            total_l += exon.end - exon.start

    for exon in exons:
        if not exon.mut_excl in ('c', 'k'):
            exon.phase = phase(total_l)
            total_l += exon.end - exon.start

    if bed.strand == "-":
        exons.sort(key = lambda ex: ex.start)


    out_str = []
    out_str.append("##gff-version 3")
    out_str.append("###")

    gene_start = min([ex.start for ex in exons]) + 1
    gene_end = max([ex.end for ex in exons])

    out_str.append("\t".join([bed.chrom,
                              ".",
                              "gene",
                              str(gene_start),
                              str(gene_end),
                              ".",
                              bed.strand,
                              ".",
                              "ID=gene:VGSC;biotype=protein_coding;Name=VGSC"]
                              )
                    )

    out_str.append("###")

    for trans in ["ck", "dl"]:
        out_str.append("\t".join([bed.chrom,
                                  ".",
                                  "mRNA",
                                  str(gene_start),
                                  str(gene_end),
                                  ".",
                                  bed.strand,
                                  ".",
                                  "ID=transcript:VGSC_" + trans + \
                                          ";Parent=gene:VGSC;Name=VGSC"\
                                           + trans + ";biotype=protein_coding"
                                 ])
                         )


        for ex_type in ["exon", "CDS"]:

            if ex_type is "exon":
                tag = "Parent=transcript:VGSC_" + trans + ";Name="
            else:
                tag = "ID=CDS:VGSC;Parent=transcript:VGSC_" + trans + ";Name="

            if trans is "dl":
                not_mut = ["c", "k"]
            else:
                not_mut = ["d", "l"]

            for exon in exons:
                if not exon.mut_excl in not_mut:
                    if ex_type is "exon":
                        phase = "."
                    else:
                        phase = str(exon.phase)
                    out_str.append("\t".join([bed.chrom,
                                              ".",
                                              ex_type,
                                              str(exon.start + 1),
                                              str(exon.end),
                                              ".",
                                              bed.strand,
                                              phase,
                                              tag + exon.name])
                         )
        out_str.append("###")

    return out_str

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Create gff3 from bed')

    parser.add_argument("bed",
                        help = "bed file")

    args = parser.parse_args()

    if not os.path.isfile(args.bed):
        print(args.bed + " does not exist!")
        sys.exit(1)

    out = create_bed_from_gff3(args.bed)

    for l in out:
        print(l)

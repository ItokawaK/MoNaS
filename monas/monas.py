#!/usr/bin/env python

import argparse
import os
import sys

try:
    from monas import genotype
    from monas import finalize_table
    from monas import make_AA_alignment
    from monas import bed2gff3
except:
    import genotype
    import finalize_table
    import make_AA_alignment
    import bed2gff3


VERSION = 2.0

def show_version(args):
    print(f'MoNaS v{VERSION} by Kentaro Itokawa')

def description(VERSION):
    return f" MoNaS (version {VERSION}) - A program genotyping VGSC genes from NGS reads."

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def main():

    parser = argparse.ArgumentParser(description=description(VERSION))
    subparsers = parser.add_subparsers(dest="command",
                         help="MoNaS command options",
                         required=True)

    parser_version = subparsers.add_parser('version', help='Show version and exit')
    parser_version.set_defaults(handler=show_version)

    # Parser run
    parser_run = subparsers.add_parser('run', help='Run all processes')
    parser_run.add_argument('-s', '--species', dest = 'species',
                        help = 'Species name. It should be same to the dirname of references.')
    parser_run.add_argument('-l', '--sample_list', dest = 'sample_list',
                        help = 'Path for a list file decribing name of sampeles and fastq files.')
    parser_run.add_argument('-t', '--max_cpu', dest = 'num_cpu',
                        type = int,
                        default = 4,
                        help = 'Maximum number of threads. [4]')
    parser_run.add_argument('-b', '--bwa_treads', dest = 'num_threads',
                        type = int,
                        default = 4,
                        help = 'Number of treads per bwa process. [4]')
    parser_run.add_argument('-o', '--out_dir', dest = 'out_dir',
                        help = 'Name of out directly. Should be new.')
    parser_run.add_argument('-r', '--ref_root', dest = 'ref_root',
                        default = SCRIPT_DIR + "/references",
                        help = 'Root directly of references. Deault = MoNaS/references')
    parser_run.add_argument('-m', '--mode', dest = 'mode',
                        default = 'ngs_dna',
                        choices = ['ngs_dna', 'ngs_rna'],
                        help = 'Analysis mode. [ngs_dna]')
    parser_run.add_argument('-c', '--variant_caller', dest = 'variant_caller',
                        default = "freebayes",
                        choices = ['freebayes', 'gatk'],
                        help = 'Variant caller to be used. Default is freebayes.')
    parser_run.add_argument('-n', '--no_clean', dest = 'do_clean',
                        action='store_false',
                        default = True,
                        help = 'Do not clean old BAM files after rmdup. Off by default.')
    parser_run.add_argument('--suppress_fullpath', dest='no_fullpath',
                        action ='store_true',
                        default= False,
                        help=argparse.SUPPRESS)
    parser_run.add_argument('--resume', dest='resume',
                        action ='store_true',
                        default= False,
                        help='Resume from existing run')
    parser_run.add_argument('-v', '--version', dest = 'show_version',
                        action='store_true',
                        default = False,
                        help = 'Show version and exit.')

    parser_run.set_defaults(handler=genotype.main)

    # Parser table
    parser_finalize = subparsers.add_parser('table',
             help='Convert a csq VCF to a table')
    parser_finalize.add_argument("out_csq",
                        help = "out_csq file path")
    parser_finalize.add_argument("ref_bed",
                        help = "ref.bed file path")
    parser_finalize.add_argument("mdom_fa",
                        help = "ref.mdom.fa file path")
    parser_finalize.set_defaults(handler=finalize_table.main)

    # Parser aln
    parser_aln = subparsers.add_parser('aln',
             help='Make AA aligment with M. domestica VGSC from DNA fasta and BED')
    parser_aln.add_argument("ref_fa",
                       help = 'Genome reference fasta')
    parser_aln.add_argument("bed",
                       help = "Reference annotation bed")
    parser_aln.add_argument("-o", "--out_aln_fasta", dest = "out_fasta_path",
                       help = "Multiple AA aligment fasta file to output.")
    parser_aln.add_argument("-t", "--translate", dest = "translate",
                       help = "Output only translation to this file.")
    parser_aln.add_argument("-m", "--mdom_path", dest = "mdom_path",
                       help = 'M. domestica aa fasta path [MoNaS/monas/misc/AAB47604.fa]',
                       default = os.path.join(SCRIPT_DIR, "misc/AAB47604.fa"))
    parser_aln.set_defaults(handler=make_AA_alignment.main)

    # Parser bed2gff
    parser_bed2gff = subparsers.add_parser('gff3',
             help='Create MoNaS compartible gff3 from BED file')
    parser_bed2gff.add_argument("bed",
                        help = "bed file")
    parser_bed2gff.set_defaults(handler=bed2gff3.main)

    args = parser.parse_args()
    args.handler(args)

if __name__ == '__main__':
    main()

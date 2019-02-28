
MoNaS
======
**An automated pipeline to genotype mosquito's voltage-gated sodium channel genes using NGS data**
### Status: <font color="Red">Version 1.0</font>

About
-------
MoNaS (**Mo**squito **Na**<sup>+</sup> channel mutation **S**earch) is an automated pipeline
conducting genotyping of voltage-gated sodium channel (VGSC) in mosquitos from NGS reads.
Basically, MoNaS is designed to utilize NGS data from genomic DNA such as targeted captured
library (SureSelect, xGen probes for instance) or RNA/cDNA such as shotgun library of PCR amplified VGSC cDNA.

Manuals
-------
Currently, we have confirmed MoNaS works in Linux OSs of both Ubuntu18 and CentOS6.
Although we have not confirmed yet, this program could also be run in Mac OS.

### Installation

MoNaS consists of several python3 scripts which do not require compilation.
However, MoNaS depends some third-party softwares described below:
- [bwa](https://github.com/lh3/bwa) v0.7.17\* (For genomic DNA data)
- [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0\* (For RNA data)
- [samtools](http://www.htslib.org/) v1.9\*
- [bcftools](http://www.htslib.org/) v1.9\*
- [freebayes](https://github.com/ekg/freebayes) v1.2.0-2\* or  [GATK 4x](https://software.broadinstitute.org/gatk/) v4.0.11.0\*
- Optionally [MUSCLE](https://www.drive5.com/muscle/) v3.8.31\*  


  \*Versions we are currently using.

Root directories of those software executables should be included in your $PATH environment variable,
or you can directly specify these paths in the **scripts/bin_path.json** file.

Additionary, you will need [biopython](https://biopython.org/) package installed in your python to run **genotype_sanger.py**, and some helper tools to create files for your own reference genomes.

### Genome references
MoNaS requires a reference genomic sequence in FASTA, and annotation for CDSs in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) files for *VGSC* gene of your species.
Accurate genome and annotation information, of course, is the most vital part of this pipeline.
Currently, MoNas includes reference files of three species of mosquitos, *Aedes aegypti* (-s Aaeg), *Aedes albopictus* (-s Aalb) and *Culex quinquefasciatus* (-s Cpip).

Each file should have name with prefix **ref** (eg. ref.ga) and be organized under
a directory as:

```
  your_species_dir/           # Arbitrary name of directory. Will be recognized as a species name (-s).
       ├- ref.fa         # Reference FASTA.
       ├- ref.fa.fai*    # FASTA index.
       ├- ref.gff3       # GFF3 annotation file for VGSC gene.
       ├- ref.bed        # BED annotation file for VGSC CDSs.
       ├- ref.dict*      # Genome dictionary used by GATK.
       ├- ref.mdom.fa    # Pair-wise aligned VGSC amino-acid sequences.
       ├- bwadb*/        # BWA indices
       │     ├- ref.amb
       │     ├- ref.ann
       │     ├- ref.bwt
       │     ├- ....
       │     
       └- hisatdb*/       # HISAT2 indices
             ├- ref.1.ht2
             ├- ref.2.ht2
             ├- .....
 *: Will be created by MoNaS automatically if absent
```
Practically, many of those files (with asterisk\*) will be created automatically by MoNaS if absent.
You will need only **re.fa** and **ref.bed** files, at least, to start analysis.
In default, MoNaS expect **your_species_dir/** locates under in the [references](https://github.com/ItokawaK/MoNaS/tree/master/references) directory of MoNaS.
You can explicitly specify the location of **your_species_dir/** by
```
-r your_reference_dir -s your_species_dir
```
in which reference files will be searched from **your_reference_dir/your_species_dir**.

- **ref.bed** annotation file.

  The 4th column (name) should be in a string with format as **ExonXX(c/d/k/l)**.
The 5th column (score) is not used.

  [example](https://raw.githubusercontent.com/ItokawaK/MoNaS/master/references/Aaeg/ref.bed)
```
3:315930000-316250000	879	1875	Exon32	0	-
3:315930000-316250000	1941	2246	Exon31	0	-
3:315930000-316250000	8708	8979	Exon30	0	-
3:315930000-316250000	9046	9292	Exon29	0	-
3:315930000-316250000	9361	9556	Exon28	0	-
3:315930000-316250000	9627	9750	Exon27	0	-
3:315930000-316250000	17891	18014	Exon26l	0	-
3:315930000-316250000	20113	20236	Exon26k	0	-
3:315930000-316250000	29958	30132	Exon25	0	-
...
3:315930000-316250000	68369	68532	Exon19d	0	-
3:315930000-316250000	69185	69348	Exon19c	0	-
3:315930000-316250000	70427	70601	Exon18	0	-
...
```

- **ref.gff3** file

  The GFF3 file interpretable by BCFtools csq.
This file can be created from **ref.bed** using [scripts/bed2gff.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/bed2gff3.py).

  [example](https://raw.githubusercontent.com/ItokawaK/MoNaS/master/references/Aaeg/ref.gff3).


- **ref.mdom.fa** file

    A multi alignment FASTA file of VGSC amino-acid (aa) sequences.
The first entry should be house fly (*M. domestica*) VGSC. The second
and third entries are aa sequences of ck and dl transcript variants
of your species, respectively.
This data will be used to convert aa position coordinate of your species to
*M. domestica* universal coordinates.
This file can also be created from **ref.fa** and **ref.bed** using [scripts/make_AA_alignment.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/make_AA_alignment.py)
(MUSCLE and biopython are required).

    [example](https://raw.githubusercontent.com/ItokawaK/MoNaS/master/references/Aaeg/ref.mdom.fa)

### Usage

```
usage:
    genotyp.py -s species_name -o out_dir_path -l list_file_path
                 (-t num_max_threads -b num_threads_per_bwa
                  -m mode[ngs_dna|ngs_rna] -r ref_dir_path
                  -v variant_caller[freebayes|gatk]
                  )


MoNaS (version 1.0) - A program genotyping VGSC genes from NGS reads.

optional arguments:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Species name. It should be same to the dirname of
                        references.
  -l SAMPLE_LIST, --sample_list SAMPLE_LIST
                        Path for a list file decribing name of sampeles and
                        fastq files.
  -t NUM_CPU, --max_cpu NUM_CPU
                        Maximum number of threads. [4]
  -b NUM_THREADS, --bwa_treads NUM_THREADS
                        Number of treads per bwa process. [4]
  -o OUT_DIR, --out_dir OUT_DIR
                        Name of out directly. Should be new.
  -r REF_ROOT, --ref_root REF_ROOT
                        Root directly of references. Deault = MoNaS/references
  -m {ngs_dna,ngs_rna}, --mode {ngs_dna,ngs_rna}
                        Analysis mode. [ngs_dna]
  -c {freebayes,gatk}, --variant_caller {freebayes,gatk}
                        Variant caller to be used. Default is freebayes.
  -n, --no_clean        Do not clean old BAM files after rmdup. Off in
                        default.
  -v, --version         Show version and exit.
```

#### Examples
1. For DNA data:
```bash
MoNaS/genotype.py  -s Aalb  -l sample_list.txt  -t 16  -o out_dir
```

2. For RNA data:
```bash
MoNaS/genotype.py  -s Aalb  -l sample_list.txt  -t 16  -o out_dir -m ngs_rna
```

#### Mandatory options

- `-s`, `--species`

  This option specifies the name of directory storing reference files of each species.
In default, this directory will be searched in the **MoNaS/references/**.
In stead, you can explicitly specify another directory using `-r, --ref_root` option (see below).

- `-l`, `--sample_list`

  Specify a path of space-delimited text file describing sample names and FASTQ paths (paried-ends or single-end) in each single row.

```example
  sample1 pe_1F.fq[.gz] pe_1R.fq[.gz]
  sample2 pe_2F.fq[.gz] pe_2R.fq[.gz]
  sample3 se_1.fq[.gz]
  sample4 pe_3F.fq[.gz] pe_3R.fq[.gz]
  ...            
```

  `sample1, sample2, ...` are arbitrary unique strings identifying each your sample. Do not include a space or characters such as /, \*, \, etc. because MoNaS will use these values for file names.
The FASTQ paths can be either absolute or relative from where you execute the `genotype.py`.

#### Other options
- `-m`, `--mode`

  Choose your sample type either from `-m ngs_dna` or `-m ngs_rna`:

     ngs_dna: NGS reads from genomic DNA (i.e. including intron) [default].
     ngs_rna: NGS reads from RNA/cDNA (i.e. introns are spliced).

- `-t`, `--max_cpu` and `-b`, `--bwa_treads`

  MoNaS uses `-t` number of threads at the same time in maximum.
  In the BWA or HISAT2 stage, several samples are processed parallelly using `-b` threads for each sample. Thus, number of BWA or HISAT2 processes running at simultaneously will be `-t` // `-b`.
  Basically, running many BWA or HISAT2 processes with small number of thread/process may
  increase speed of analysis but increases memory usage especially if you use
  large genome reference (e.g. full genome).

- `-r`, `--ref_root`

  You can specify, path of root directory where the species directory will be searched.

- `-c`, `--variant_caller`

  You can select variant caller program to be used (freebayes or gatk).
  Currently, we have tested MoNaS extensively with freebayes more than gatk.  

Pipeline detail
--------------

1. MoNaS maps NGS reads to reference genome of each mosquito species with `bwa mem` for DNA data or `hisat2` for RNA data.

1. The resulted bam files are sorted, removed PCR duplicates, and indexed with `samtools sort`, `rmdup` and `index`, respectively.

1. Each indexed bam file are processed with `freebayes` or `gatk HaplotypeCaller`. **ref.bed** will be used to restrict regions to be analyzed. `freebayes` processes all bam files as single run, but multiprocessed by dividing the region of interest into many sub-regions (exons) which will be integrated in single **out.vcf** file at the end. `gatk HaplotypeCaller`, on the other hand, processes only single bam at once, creating vcf files for each sample.

1. Annotates the **out.vcf** for amino acid changes by `bcftools csq -p a -l` using information in **ref.gff3** resulting in **out_csq.vcf**.

1. Finally, human-friendly table **table_with_Mdomcoord.tsv** describing amino acid changes and its corresponding AA positions in *Musca domestica* will be generated from **out_csq.vcf** files.

Output
------
The output directry will look like:

```bash
  out_dir/ # your specified name
       ├- BAMs/     # Sorted bam files. Will be removed after samtools rmdupped in default.
       │    ├- sample1.bam
       │    ├- sample2.bam
       │    ├- ...
       │
       ├- BAMs_rmdup/　# Indexed bam files after remove PCR duplicates.
       │    ├- sample1.bam
       │    ├- sample1.bam.bai
       │    ├- sample2.bam
       │    ├- ...
       │       
       ├- VCFs/ # indexed vcf files for each individual by gatk
       │     ├- sample1.vcf
       │     ├- sample1.vcf.idx
       │     ├- sample2.vcf
       │     ├- ....
       │     
       ├- out.vcf # vcf file for multiple samples by freebayes
       ├- out_csq.vcf # vcf file for multiple samples with csq tag
       └- table_with_Mdomcoord.tsv # list of mutations and AA changes with M. domestica AA number



```

table_with_Mdomcoord.tsv
------
The final output table, **table_with_Mdomcoord.tsv**, will look like as below.

```
#ID     CHROM   POS     REF_ALLELE      ALT_ALLELE      QUAL    GENOTYPE        AA_CHANGE       AA_CHANGE_MDOM  KDR_EVIDENCE    AD      EXON
Aalb-Yona-02    MNAF02001058.1  1930245 T       C       6445.79 T/C     wild/14S        wild/synonymous NA/NA   287,249 Exon1
Aalb-Okayama-02 MNAF02001058.1  1978801 G       T       10826.2 G/T     wild/74P        wild/synonymous NA/NA   163,172 Exon3
Aalb-Yona-04    MNAF02001058.1  1978801 G       T       10826.2 G/T     wild/74P        wild/synonymous NA/NA   263,239 Exon3
Aalb-Okayama-02 MNAF02001058.1  1978810 A       G       10877.2 A/G     wild/77T        wild/synonymous NA/NA   163,173 Exon3
...
Aalb-Yona-08    MNAF02001058.1  2179316 CATC    TATT    507016.0        TATT/TATT       1335VI/1335VI   synonymous/synonymous   NA/NA   0,286   Exon24
Aalb-SP-02      MNAF02001058.1  2179316 CATC    TATT    507016.0        TATT/TATT       1335VI/1335VI   synonymous/synonymous   NA/NA   0,1307  Exon24
Aalb-SP-04      MNAF02001058.1  2179316 CATC    TATT    507016.0        TATT/TATT       1335VI/1335VI   synonymous/synonymous   NA/NA   0,1008  Exon24
Aalb-SP-07      MNAF02001058.1  2179316 CATC    TATT    507016.0        TATT/TATT       1335VI/1335VI   synonymous/synonymous   NA/NA   0,1272  Exon24
Aalb-SP-06      MNAF02001058.1  2179316 CATC    TATT    507016.0        TATT/TATT       1335VI/1335VI   synonymous/synonymous   NA/NA   0,1104  Exon24
...
Aalb-SP-03      MNAF02001058.1  2207911 T       G       231649.0        G/G     F1574C/F1574C   F1534C/F1534C   Strong/Strong   0,906   Exon29
Aalb-SP-01      MNAF02001058.1  2207911 T       G       231649.0        G/G     F1574C/F1574C   F1534C/F1534C   Strong/Strong   1,934   Exon29
Aalb-SP-02      MNAF02001058.1  2207911 T       G       231649.0        G/G     F1574C/F1574C   F1534C/F1534C   Strong/Strong   1,990   Exon29
```
Column 1: Sample ID

Column 2: Name of scaffold

Column 3: Nucleotide position of variant

Column 4: Reference allele at this site

Column 5: One or set of alternative alleles at this site

Column 6: Quality score from variant caller

Column 7: Genotype of this sample

Column 8: Annotated AA change **in reference AA position**

Column 9: Annotated AA change **in *M. domestica* AA position**.

Collum 10: Whether AA change is known *kdr* mutation or not.

  - AA variant which is confirmed to decrease VGSC sensitivity to pyrethroid
  in electrophysiological assay will be described "Strong".
  Those with weaker evidences are described as "Weak".
  Other unknown AA variants will be described as "Unknown".
  Known kdr AA varoamts are listed in [scripts/kdr_list.json](https://github.com/ItokawaK/MoNaS/blob/master/scripts/kdr_list.json) which will be updated as possible as regularly.

  **Note that the list would not always be complete!**.

Column 11: Read depth for each allele

Column 12: Exon where this variant belongs

### Aanlyzing Sanger sequence reads
  Although it does not seem the best approach, reads from Sanger sequence technology could be analyzed in
MoNaS pipeline by regarding those Sanger reads as NGS reads. [genotype_sanger.py](https://github.com/ItokawaK/MoNaS/blob/master/genotype_sanger.py) is a wrapper script to conduct
chopping input Sanger reads into 150 bp short reads (~ x5 coverage) with fake quality values, writing fastq and sample list files,
and then executing MoNaS for those data.

  As any sequecing errors existing in input Sanger reads will be considered as **true variants**, it is important to trim low-quality regions in advance. Also, ambiguous nucleotide codes (R, Y, S, etc...) are not supported yet (ToDo).

```bash
MoNaS/genotype_sanger.py -s Aalb -t 16 -o out_table.tsv sanger_reads.fa
```


Helper tools
------

MoNaS includes some tools assisting creation of new reference annotation file for species of your interest.

- [scripts/bed2gff.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/bed2gff3.py)

  This script creates a gff3 file which is interpretable by bcftools csq from a bed file describing *VGSC* CDSs.

```
usage: bed2gff3.py [-h] bed

Create gff3 from bed

positional arguments:
  bed         bed file

optional arguments:
  -h, --help  show this help message and exit
```

- [scripts/make_AA_alignment.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/make_AA_alignment.py)

  This script translates genome to VGSC protein using CDS information described in bed file, then conducts
pairwise alignment with VGSC in *M. domestica* which is usable as **ref.mdom.fa**.
The script also reports mismatched AA between your reference VGSC and *M. domestica* in **VGSC  OUT_FASTA_PATH.info** with notification for potential kdr(s) listed in [scripts/kdr_list.json](https://github.com/ItokawaK/MoNaS/blob/master/scripts/kdr_list.json) if found.

```
usage: make_AA_alignment.py [-h] [-o OUT_FASTA_PATH] [-t TRANSLATE]
                            [-m MDOM_PATH]
                            ref_fa bed

Make AA aligment with M. domestica VGSCfrom dna fasta and bed

positional arguments:
  ref_fa                Reference fasta
  bed                   Reference annotation bed

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_FASTA_PATH, --out_aln_fasta OUT_FASTA_PATH
                        Multiple AA aligment fasta file to output.
  -t TRANSLATE, --translate TRANSLATE
                        Output only translation to this file.
  -m MDOM_PATH, --mdom_path MDOM_PATH
                        M. domestica aa fasta path [MoNaS/misc/AAB47604.fa]
```

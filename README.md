
MoNaS 
======
**An automated pipeline to genotype mosquito's voltage-gated sodum channel genes using NGS data**
### Status: <font color="Red">Under construction. Sorry!</font>

MoNaS (**Mo**squito **Na**<sup>+</sup> channel mutation **S**earch) is an automated pipeline to assist genotyping
voltage-gated sodium channel (VGSC) in mosquitos.
Basically, MoNaS is designed to handle NGS data from genomic DNA such as targeted captured library (SureSelect, xGen probes for instance) or RNA/cDNA such as shotgun library of PCR amplified VGSC cDNA. 

Manuals
-------
Currently, we have confirmed MoNaS works in both Ubuntu18 and CentOS6. It is possible to run in Mac, perhaps. 

### Installation

MoNaS consists of several python3 scripts which do not require compilation.
However, MoNaS requires some third party softwares: 
- [bwa](https://github.com/lh3/bwa) v0.7.17 (For genomic DNA data)
- [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0 (For RNA data)
- [samtools](http://www.htslib.org/) v1.9
- [bcftools](http://www.htslib.org/) v1.9
- [GATK 4x](https://software.broadinstitute.org/gatk/) v4.0.11.0
 \*Versions we used. 

Root directories of those sofwares should be included in your $PATH, or you can directly
specify them in the **scripts/bin_path.json** file.

### Genome references
MoNaS requires a reference fasta sequence (ref.fa), and annotation files (ref.gff3, re.bed) for each mosquito species.
Actually, accurate information for VGSC exons (in gff3 format) is the heart of this pipeline. Currently, we have
prepared reference files of three species of mosquitos, *Aedes aegypti*, *Aedes albopictus* and *Culex quinquefasciatus* which
will be uploaded somewhere (please wait for a while).

### Usage

1. For DNA data: 
```bash
MoNaS/genotype.py  -s species_name  -l sample_list.txt  -t num_cpu  -o out_dir
```

2. For RNA data:
```bash
MoNaS/genotype.py  -s species_name  -l sample_list.txt  -t num_cpu  -o out_dir -m ngs_rna
```
### Option details

- `-s`, `--species`
  
This option specifies the name of directory storing reference files of each species. Each file should have file names with
prefix ref (eg. ref.ga) and be organiazed within this directory as:

```bash
  species_dir/ # arbitrary name
       ├- ref.fa      # reference fasta
       ├- ref.fa.fai  # fasta index
       ├- ref.gff3    # gff3 annotation for VGSC CDSs
       ├- ref.bed     # bed annotation for VGSC CDSs
       ├- ref.dict    # genome dictionary for GATK
       ├- ref.mdom.fa # pair-wise aligned VGSC AA seqs of mosquito and house-fly (genbank_id: AAB47604)
       ├- bwadb/ # bwa indicies
       │     ├- ref.amb
       │     ├- ref.ann
       │     ├- ref.bwt
       │     ├- ....
       │     
       └-hisatdb/ # hisat2 indicies
             ├-ref.1.ht2
             ├-ref.2.ht2
             ├- .....
             
```
In deafault, the **species_dir/** will be searched in the **MoNaS/references/**. In stead, 
you can explicitly specify another directory using `-r, --ref_root` option. If the program could not find
**bwadb/**, **hisatdb/**, **ref.fa.fai** or **ref.dict** within the reference directory, `getnotype.py` will automaticaly 
try to create them.

- `-l`, `--sample_list`

 Specify a path of space-delimited text file describing sample names and FASTQ paths (paried-ends or single-end) in each single row.

```example
  sample1 pe_1F.fq[.gz] pe_1R.fq[.gz]
  sample2 pe_2F.fq[.gz] pe_2R.fq[.gz]
  sample3 se_1.fq[.gz]
  sample4 pe_3F.fq[.gz] pe_3R.fq[.gz]
  ...            
```

`sample1, sample2, ...` are arbitrary unique strings identifying your each sample. Do not include a space or characters such as /, \*, \, etc. because MoNaS will use these values for file names. The FASTQ paths can be relative from where you call `genotype.py`.

- `-m`, `--mode`

Choose your sample type either from `-m ngs_dna` or `-m ngs_rna`:

     ngs_dna: NGS reads from genomic DNA (i.e. including intron) [default].
     ngs_rna: NGS reads from RNA/cDNA (i.e. introns are spliced).

- fasta



Pipeline detail
--------------

1. MoNaS maps NGS reads to reference genome of each mosquito species with `bwa mem` for DNA data or `hisat2` for RNA data.

1. The resulted bam files are sorted, removed PCR duplicates, and indexed with `samtools sort`, `rmdup` and `index`, respectively. 

1. Each indexed bam file are processed with `gatk HaplotypeCaller`. **ref.bed** will be used to restrict regions analyzed.  

1. `bcftools csq` annotate amino acid changes using information in **ref.gff3**.

1. Finally, hunam-friendly table decribing amino acid changes and its corresponding AA positions in *Musca domestica* will be
generated from vcf files.

Output
------
The output data will look like:

```bash
  out_dir/ # your specified name
       ├- BAMs/     # indexed and sorted bam files
       │    ├- sample1.bam
       │    ├- sample1.bam.bai
       │    ├- sample2.bam
       │    ├- ...
       │
       ├- VCFs/ # indexed vcf files
       │     ├- sample1.vcf
       │     ├- sample1.vcf.idx
       │     ├- sample2.vcf
       │     ├- ....
       │     
       ├- out_table # 
       └- table_with_Mdomcoord.tsv # list of mutations and AA changes with M. domestica AA number
            
```



MoNaS 
======
**A pipeline to genotype mosquito's voltage-gated sodum channel genes using NGS data**
### Status: <font color="Red">Under construction. Sorry!</font>

MoNaS (**Mo**squito **Na**<sup>+</sup> channel mutation **S**earch) is a pipelin tool to genotype
voltage-gated sodium channel (VGSC) in mosquitos.
 Basically, MoNaS is designed to
handle NGS data from targeted captured library (SureSelect, xGen probes for instance), but we are planning to let this
tool to treat data of Sanger sequence and RNA-seq.

Manuals
-------

### Instlation

Currently, we have confirmed MoNaS works in Ubuntu18 and CentOS6. It may also work in Mac. 

MoNaS consists of several python3 scripts which do not require compilation.
However, MoNaS requires some third party softwares described below.
- [bwa](https://github.com/lh3/bwa) (For genomic DNA data)
- [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) (For RNA data)
- [samtools](http://www.htslib.org/)
- [bcftools](http://www.htslib.org/)
- [GATK 4x](https://software.broadinstitute.org/gatk/)

Root directories of those sofwares should be included in your $PATH, or you can directly
specify them in the **scripts/bin_path.json** file.

MoNaS also requires a reference fasta sequence (ref.fa), and annotation files (ref.gff3, re.bed) 
for each mosquito species. Currently, reference of MoNaS includes three species of mosquitos, *Aedes aegypti*, 
*Aedes albopictus* and *Culex quinquefasciatus*. 

### Usage

1. For DNA data: 
```bash
MoNaS/genotype.py  -s species_name  -l sample_list.txt  -t num_cpu  -o out_dir
```

2. For RNA data:
```bash
MoNaS/genotype.py  -s species_name  -l sample_list.txt  -t num_cpu  -o out_dir -m ngs_rna
```

- `-s`, `--species`
  
The name directory storing reference files of each species. Each file should have file names with prefix ref (eg. ref.ga) and be
organiazed under this directory as:

```bash
  species_name/ (arbitrary name)
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
In deafault, the reference directory will be searched for the **references/** directory in the script directory. In stead, 
you can specify another directory explicitly with `-r, --ref_root` option. If the program could not find
**bwadb/**, **hisatdb/**, **ref.fa.fai** or **ref.dict** within the reference directory, `getnotype.py` will automaticaly 
try to create them.


- `-l`, `--sample_list`

  A space delimited text file describing sample names and FASTQ paths (paried-ends or single-end) in single row.

```
  sample1 pe_1F.fq[.gz] pe_1R.fq[.gz]
  sample2 pe_2F.fq[.gz] pe_2R.fq[.gz]
  sample3 se_1.fq[.gz]
  sample4 pe_3F.fq[.gz] pe_3R.fq[.gz]
  ...            
```

The *sample_names* are arbitrary unique strings identifying your each sample. Do not include a space or characters such as / * \% etc.., because MoNaS will use this value as a file name. The FASTQ paths can be relative from where you call `genotype.py`.

- `-m`, `--mode`

Choose your sample type among:

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


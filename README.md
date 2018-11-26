
# MoNaS - <u>Mo</u>squito <u>Na</u><sup>+</sup> channel mutation <u>S</u>earch
==============

### Status: <font color="Red">Under construction. Sorry! </font>

MoNaS is a pipelin tool to genotype voltage-gated sodium channel (VGSC) in mosquitos.
Currently reference of MoNaS includes three species of mosquitos, *Aedes aegypti*, 
*Aedes albopictus* and *Culex quinquefasciatus*. Basically, MoNaS is designed to
handle NGS data from targeted captured library (SureSelect, xGen probes for instance), but we are planning to let this
tool to treat data from sanger sequence or RNA-seq.



Manuals
-------

### Instlation

Currently, we have confirmed MoNaS works in Ubuntu18 and CentOS6. It may work also
in Mac. 

MoNaS consists of several python3 scripts which do not require compilation.
However, MoNaS requires some third party softwares described below are
installed.
- bwa
- samtools
- bcftools
- GATK

Root directories of those sofwares should be included in your $PATH, or you can directly
describe them in the **configuration.py** script file.

MoNaS also requires a reference genome sequence (fasta), its bwa indices and annotation files 
(gff3, bed) for each mosquito species. The fasta file should be pre-indexed by *samtools faidx*. 

### Usage

```bash
[MoNaS_v1.* directory path]/genotype.py  -s ref_root_dir_name  -l sample_list  -t num_cpu  -o out_dir
```
-r species_ref_root_dir
  
The directory storing reference files of each species. These files should have file names with prefix ref (eg. ref.ga) and be
organiazed as:

```bash
  ref_dir_name (arbitrary name)
       ├- ref.fa      # reference fasta
       ├- ref.fa.fai  # fasta index
       ├- ref.gff3    # gff3 annotation for VGSC CDSs
       ├- ref.bed     # bed annotation for VGSC CDSs
       └- bwadb/
            ├- ref.amb
            ├- ref.ann
            ├- ref.bwt
            ├- ref.pac
            └- ref.sa # bwa indices for ref.fa
       
```
In deafault, the reference directory is expexted to locate in the "references/" directory in the script directory. In stead, 
you can explicitly assign the root directory of references with the -r, --ref_root option.


- sample_list

  A space delimited text file describing sample names and FASTQ of paried-ends or single-end .

```
  sample_name_1 fastq_1F_path.fq[.gz] fastq_1R_path.fq[.gz]
  sample_name_2 fastq_2F_path.fq[.gz] fastq_2R_path.fq[.gz]
  sample_name_3 fastq_3_path.fq[.gz]
  sample_name_4 fastq_4F_path.fq[.gz] fastq_4R_path.fq[.gz]
```

The *sample_names* are arbitrary strings but do not include characters such as / * \% etc.., and they should be unique because MoNaS will use *sample_names* as file names. The FASTQ path can be relative from where you run MoNaS.

Pipeline detail
--------------

MoNaS first maps NGS reads to reference genome of each mosquito species with *bwa mem*. 
Then, the resulted bam file is sorted, removed PCR duplicates, and indexed with *samtools* 
functions. The indexed bam file is processed with *GATK* for variant calling for the VGSC gene,
and then, amino acid changes are annotated with *bcftools csq*. Finally, hunam-friendly table
decribing amino acid changes and its corresponding AA positions in *Musca domestica* will be
generated.


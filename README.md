
MoNaS 
======
*A pipeline to genotype voltage-gated sodum channel in mosquitos*
### Status: <font color="Red">Under construction. Sorry!</font>

MoNaS (**Mo**squito **Na**<sup>+</sup> channel mutation **S**earch) is a pipelin tool to genotype
voltage-gated sodium channel (VGSC) in mosquitos.
 Basically, MoNaS is designed to
handle NGS data from targeted captured library (SureSelect, xGen probes for instance), but we are planning to let this
tool to treat data of Sanger sequence and RNA-seq.

Manuals
-------

### Instlation

Currently, we have confirmed MoNaS works in Ubuntu18 and CentOS6. It may work also
in Mac. 

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

```bash
MoNaS/genotype.py  -s species_name  -l sample_list.txt  -t num_cpu  -o out_dir
```
- `-s, --species`
  
The name directory storing reference files of each species. Each file should have file names with prefix ref (eg. ref.ga) and be
organiazed under this directory as:

```bash
  species_name/ (arbitrary name)
       ├- ref.fa      # reference fasta
       ├- ref.fa.fai  # fasta index
       ├- ref.gff3    # gff3 annotation for VGSC CDSs
       ├- ref.bed     # bed annotation for VGSC CDSs
       ├- ref.dict    # genome dictionary for GATK
       ├- bwadb/ # bwa indicies
       |     ├- ref.amb
       |     ├- ref.ann
       |     ├- ref.bwt
       |     ├- ....
       |     
       └-hisatdb/ # hisat2 indicies
             ├-ref.1.ht2
             ├-ref.2.ht2
             ├- .....
             
```
In deafault, the reference directory will be searched for the **references/** directory in the script directory. In stead, 
you can explicitly specify another directory of references with the -r, --ref_root option. If the program could not find
**bwadb/**, **hisatdb/**, **ref.fa.fai** or **ref.dict**, `getnotype.py` will automaticaly create them.


- `-l, --sample_list`

  A space delimited text file describing sample names and specifying FASTQ paths of paried-ends or single-end .

```
  sample_name_1 pe_1F.fq[.gz] pe_1R.fq[.gz]
  sample_name_2 pe_2F.fq[.gz] pe_2R.fq[.gz]
  sample_name_3 se_1.fq[.gz]
  sample_name_4 pe_3F.fq[.gz] pe_3R.fq[.gz]
  ...            
```

The *sample_names* are arbitrary unique strings identifying your each sample. Do not include a space or characters such as / * \% etc.., because MoNaS will use this value as a file name. The FASTQ path can be relative from where you run MoNaS.

- `-m, --mode`

Choose your sample type among:

 + `ngs_dna`: NGS reads from genomic DNA (i.e. including intron) \[defaul\].
 
 + `ngs_rna`: NGS reads from RNA/cDNA (i.e. introns are spliced).

- fasta



Pipeline detail
--------------

MoNaS first maps NGS reads to reference genome of each mosquito species with *bwa mem*. 
Then, the resulted bam file is sorted, removed PCR duplicates, and indexed with *samtools* 
functions. The indexed bam file is processed with *GATK* for variant calling for the VGSC gene,
and then, amino acid changes are annotated with *bcftools csq*. Finally, hunam-friendly table
decribing amino acid changes and its corresponding AA positions in *Musca domestica* will be
generated.


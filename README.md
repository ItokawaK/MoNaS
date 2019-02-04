
MoNaS
======
**An automated pipeline to genotype mosquito's voltage-gated sodium channel genes using NGS data**
### Status: <font color="Red">Under construction. Sorry!</font>

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
Currently, we have reference files of three species of mosquitos, *Aedes aegypti*, *Aedes albopictus* and *Culex quinquefasciatus* which will be uploaded somewhere (please wait for a while).

Each file should have name with prefix **ref** (eg. ref.ga) and be organized under
a directory as

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
In default, MoNaS expect **your_species_dir/** locates under in the **references/** directory of MoNaS.
You can explicitly specify the location of **your_species_dir/** by
```
-r your_reference_dir -s your_species_dir
```
in which reference files will be searched from **your_reference_dir/your_species_dir**.

- Example of **ref.bed** annotation file.

4-th column should be in format **ExonXX(c/d/k/l)**.

```
3:315930000-316250000	879	1875	Exon32
3:315930000-316250000	1941	2246	Exon31
3:315930000-316250000	8708	8979	Exon30
3:315930000-316250000	9046	9292	Exon29
3:315930000-316250000	9361	9556	Exon28
3:315930000-316250000	9627	9750	Exon27
3:315930000-316250000	17891	18014	Exon26l
3:315930000-316250000	20113	20236	Exon26k
3:315930000-316250000	29958	30132	Exon25
3:315930000-316250000	30224	30490	Exon24
3:315930000-316250000	30555	30775	Exon23
3:315930000-316250000	42191	42403	Exon22
3:315930000-316250000	53517	53764	Exon21
3:315930000-316250000	53997	54185	Exon20
3:315930000-316250000	68369	68532	Exon19d
3:315930000-316250000	69185	69348	Exon19c
3:315930000-316250000	70427	70601	Exon18
3:315930000-316250000	70707	70985	Exon17
3:315930000-316250000	84047	84086	Exon16.5
3:315930000-316250000	84426	84495	Exon16
3:315930000-316250000	84558	84658	Exon15
3:315930000-316250000	99303	99506	Exon14
3:315930000-316250000	100295	100479	Exon13
3:315930000-316250000	100610	100673	Exon12
3:315930000-316250000	137698	137974	Exon11
3:315930000-316250000	150630	150795	Exon10
3:315930000-316250000	151034	151179	Exon9
3:315930000-316250000	152010	152071	Exon8
3:315930000-316250000	152157	152370	Exon7
3:315930000-316250000	152469	152561	Exon6
3:315930000-316250000	171848	171977	Exon5
3:315930000-316250000	199583	199789	Exon4
3:315930000-316250000	240354	240510	Exon3
3:315930000-316250000	273651	273684	Exon2
3:315930000-316250000	313935	314082	Exon1
```

- Example of **ref.gff3** file
The GFF3 file interpretable by BCFtools csq. This file can be created from **ref.bed**
using [scripts/bed2gff.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/bed2gff3.py).

```
##gff-version 3
###
3:315930000-316250000	.	gene	880	314082	.	-	.	ID=gene:VGSC;biotype=protein_coding;Name=VGSC
###
3:315930000-316250000	.	mRNA	880	314082	.	-	.	ID=transcript:VGSC_ck;Parent=gene:VGSC;Name=VGSCck;biotype=protein_coding
3:315930000-316250000	.	exon	880	1875	.	-	.	Parent=transcript:VGSC_ck;Name=Exon32
3:315930000-316250000	.	exon	1942	2246	.	-	.	Parent=transcript:VGSC_ck;Name=Exon31
3:315930000-316250000	.	exon	8709	8979	.	-	.	Parent=transcript:VGSC_ck;Name=Exon30
3:315930000-316250000	.	exon	9047	9292	.	-	.	Parent=transcript:VGSC_ck;Name=Exon29
3:315930000-316250000	.	exon	9362	9556	.	-	.	Parent=transcript:VGSC_ck;Name=Exon28
3:315930000-316250000	.	exon	9628	9750	.	-	.	Parent=transcript:VGSC_ck;Name=Exon27
3:315930000-316250000	.	exon	20114	20236	.	-	.	Parent=transcript:VGSC_ck;Name=Exon26k
3:315930000-316250000	.	exon	29959	30132	.	-	.	Parent=transcript:VGSC_ck;Name=Exon25
3:315930000-316250000	.	exon	30225	30490	.	-	.	Parent=transcript:VGSC_ck;Name=Exon24
3:315930000-316250000	.	exon	30556	30775	.	-	.	Parent=transcript:VGSC_ck;Name=Exon23
3:315930000-316250000	.	exon	42192	42403	.	-	.	Parent=transcript:VGSC_ck;Name=Exon22
3:315930000-316250000	.	exon	53518	53764	.	-	.	Parent=transcript:VGSC_ck;Name=Exon21
3:315930000-316250000	.	exon	53998	54185	.	-	.	Parent=transcript:VGSC_ck;Name=Exon20
3:315930000-316250000	.	exon	69186	69348	.	-	.	Parent=transcript:VGSC_ck;Name=Exon19c
3:315930000-316250000	.	exon	70428	70601	.	-	.	Parent=transcript:VGSC_ck;Name=Exon18
3:315930000-316250000	.	exon	70708	70985	.	-	.	Parent=transcript:VGSC_ck;Name=Exon17
3:315930000-316250000	.	exon	84048	84086	.	-	.	Parent=transcript:VGSC_ck;Name=Exon16.5
3:315930000-316250000	.	exon	84427	84495	.	-	.	Parent=transcript:VGSC_ck;Name=Exon16
3:315930000-316250000	.	exon	84559	84658	.	-	.	Parent=transcript:VGSC_ck;Name=Exon15
3:315930000-316250000	.	exon	99304	99506	.	-	.	Parent=transcript:VGSC_ck;Name=Exon14
3:315930000-316250000	.	exon	100296	100479	.	-	.	Parent=transcript:VGSC_ck;Name=Exon13
3:315930000-316250000	.	exon	100611	100673	.	-	.	Parent=transcript:VGSC_ck;Name=Exon12
3:315930000-316250000	.	exon	137699	137974	.	-	.	Parent=transcript:VGSC_ck;Name=Exon11
3:315930000-316250000	.	exon	150631	150795	.	-	.	Parent=transcript:VGSC_ck;Name=Exon10
3:315930000-316250000	.	exon	151035	151179	.	-	.	Parent=transcript:VGSC_ck;Name=Exon9
3:315930000-316250000	.	exon	152011	152071	.	-	.	Parent=transcript:VGSC_ck;Name=Exon8
3:315930000-316250000	.	exon	152158	152370	.	-	.	Parent=transcript:VGSC_ck;Name=Exon7
3:315930000-316250000	.	exon	152470	152561	.	-	.	Parent=transcript:VGSC_ck;Name=Exon6
3:315930000-316250000	.	exon	171849	171977	.	-	.	Parent=transcript:VGSC_ck;Name=Exon5
3:315930000-316250000	.	exon	199584	199789	.	-	.	Parent=transcript:VGSC_ck;Name=Exon4
3:315930000-316250000	.	exon	240355	240510	.	-	.	Parent=transcript:VGSC_ck;Name=Exon3
3:315930000-316250000	.	exon	273652	273684	.	-	.	Parent=transcript:VGSC_ck;Name=Exon2
3:315930000-316250000	.	exon	313936	314082	.	-	.	Parent=transcript:VGSC_ck;Name=Exon1
3:315930000-316250000	.	CDS	880	1875	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon32
3:315930000-316250000	.	CDS	1942	2246	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon31
3:315930000-316250000	.	CDS	8709	8979	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon30
3:315930000-316250000	.	CDS	9047	9292	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon29
3:315930000-316250000	.	CDS	9362	9556	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon28
3:315930000-316250000	.	CDS	9628	9750	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon27
3:315930000-316250000	.	CDS	20114	20236	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon26k
3:315930000-316250000	.	CDS	29959	30132	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon25
3:315930000-316250000	.	CDS	30225	30490	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon24
3:315930000-316250000	.	CDS	30556	30775	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon23
3:315930000-316250000	.	CDS	42192	42403	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon22
3:315930000-316250000	.	CDS	53518	53764	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon21
3:315930000-316250000	.	CDS	53998	54185	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon20
3:315930000-316250000	.	CDS	69186	69348	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon19c
3:315930000-316250000	.	CDS	70428	70601	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon18
3:315930000-316250000	.	CDS	70708	70985	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon17
3:315930000-316250000	.	CDS	84048	84086	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon16.5
3:315930000-316250000	.	CDS	84427	84495	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon16
3:315930000-316250000	.	CDS	84559	84658	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon15
3:315930000-316250000	.	CDS	99304	99506	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon14
3:315930000-316250000	.	CDS	100296	100479	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon13
3:315930000-316250000	.	CDS	100611	100673	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon12
3:315930000-316250000	.	CDS	137699	137974	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon11
3:315930000-316250000	.	CDS	150631	150795	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon10
3:315930000-316250000	.	CDS	151035	151179	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon9
3:315930000-316250000	.	CDS	152011	152071	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon8
3:315930000-316250000	.	CDS	152158	152370	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon7
3:315930000-316250000	.	CDS	152470	152561	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon6
3:315930000-316250000	.	CDS	171849	171977	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon5
3:315930000-316250000	.	CDS	199584	199789	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon4
3:315930000-316250000	.	CDS	240355	240510	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon3
3:315930000-316250000	.	CDS	273652	273684	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon2
3:315930000-316250000	.	CDS	313936	314082	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_ck;Name=Exon1
###
3:315930000-316250000	.	mRNA	880	314082	.	-	.	ID=transcript:VGSC_dl;Parent=gene:VGSC;Name=VGSCdl;biotype=protein_coding
3:315930000-316250000	.	exon	880	1875	.	-	.	Parent=transcript:VGSC_dl;Name=Exon32
3:315930000-316250000	.	exon	1942	2246	.	-	.	Parent=transcript:VGSC_dl;Name=Exon31
3:315930000-316250000	.	exon	8709	8979	.	-	.	Parent=transcript:VGSC_dl;Name=Exon30
3:315930000-316250000	.	exon	9047	9292	.	-	.	Parent=transcript:VGSC_dl;Name=Exon29
3:315930000-316250000	.	exon	9362	9556	.	-	.	Parent=transcript:VGSC_dl;Name=Exon28
3:315930000-316250000	.	exon	9628	9750	.	-	.	Parent=transcript:VGSC_dl;Name=Exon27
3:315930000-316250000	.	exon	17892	18014	.	-	.	Parent=transcript:VGSC_dl;Name=Exon26l
3:315930000-316250000	.	exon	29959	30132	.	-	.	Parent=transcript:VGSC_dl;Name=Exon25
3:315930000-316250000	.	exon	30225	30490	.	-	.	Parent=transcript:VGSC_dl;Name=Exon24
3:315930000-316250000	.	exon	30556	30775	.	-	.	Parent=transcript:VGSC_dl;Name=Exon23
3:315930000-316250000	.	exon	42192	42403	.	-	.	Parent=transcript:VGSC_dl;Name=Exon22
3:315930000-316250000	.	exon	53518	53764	.	-	.	Parent=transcript:VGSC_dl;Name=Exon21
3:315930000-316250000	.	exon	53998	54185	.	-	.	Parent=transcript:VGSC_dl;Name=Exon20
3:315930000-316250000	.	exon	68370	68532	.	-	.	Parent=transcript:VGSC_dl;Name=Exon19d
3:315930000-316250000	.	exon	70428	70601	.	-	.	Parent=transcript:VGSC_dl;Name=Exon18
3:315930000-316250000	.	exon	70708	70985	.	-	.	Parent=transcript:VGSC_dl;Name=Exon17
3:315930000-316250000	.	exon	84048	84086	.	-	.	Parent=transcript:VGSC_dl;Name=Exon16.5
3:315930000-316250000	.	exon	84427	84495	.	-	.	Parent=transcript:VGSC_dl;Name=Exon16
3:315930000-316250000	.	exon	84559	84658	.	-	.	Parent=transcript:VGSC_dl;Name=Exon15
3:315930000-316250000	.	exon	99304	99506	.	-	.	Parent=transcript:VGSC_dl;Name=Exon14
3:315930000-316250000	.	exon	100296	100479	.	-	.	Parent=transcript:VGSC_dl;Name=Exon13
3:315930000-316250000	.	exon	100611	100673	.	-	.	Parent=transcript:VGSC_dl;Name=Exon12
3:315930000-316250000	.	exon	137699	137974	.	-	.	Parent=transcript:VGSC_dl;Name=Exon11
3:315930000-316250000	.	exon	150631	150795	.	-	.	Parent=transcript:VGSC_dl;Name=Exon10
3:315930000-316250000	.	exon	151035	151179	.	-	.	Parent=transcript:VGSC_dl;Name=Exon9
3:315930000-316250000	.	exon	152011	152071	.	-	.	Parent=transcript:VGSC_dl;Name=Exon8
3:315930000-316250000	.	exon	152158	152370	.	-	.	Parent=transcript:VGSC_dl;Name=Exon7
3:315930000-316250000	.	exon	152470	152561	.	-	.	Parent=transcript:VGSC_dl;Name=Exon6
3:315930000-316250000	.	exon	171849	171977	.	-	.	Parent=transcript:VGSC_dl;Name=Exon5
3:315930000-316250000	.	exon	199584	199789	.	-	.	Parent=transcript:VGSC_dl;Name=Exon4
3:315930000-316250000	.	exon	240355	240510	.	-	.	Parent=transcript:VGSC_dl;Name=Exon3
3:315930000-316250000	.	exon	273652	273684	.	-	.	Parent=transcript:VGSC_dl;Name=Exon2
3:315930000-316250000	.	exon	313936	314082	.	-	.	Parent=transcript:VGSC_dl;Name=Exon1
3:315930000-316250000	.	CDS	880	1875	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon32
3:315930000-316250000	.	CDS	1942	2246	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon31
3:315930000-316250000	.	CDS	8709	8979	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon30
3:315930000-316250000	.	CDS	9047	9292	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon29
3:315930000-316250000	.	CDS	9362	9556	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon28
3:315930000-316250000	.	CDS	9628	9750	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon27
3:315930000-316250000	.	CDS	17892	18014	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon26l
3:315930000-316250000	.	CDS	29959	30132	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon25
3:315930000-316250000	.	CDS	30225	30490	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon24
3:315930000-316250000	.	CDS	30556	30775	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon23
3:315930000-316250000	.	CDS	42192	42403	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon22
3:315930000-316250000	.	CDS	53518	53764	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon21
3:315930000-316250000	.	CDS	53998	54185	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon20
3:315930000-316250000	.	CDS	68370	68532	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon19d
3:315930000-316250000	.	CDS	70428	70601	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon18
3:315930000-316250000	.	CDS	70708	70985	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon17
3:315930000-316250000	.	CDS	84048	84086	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon16.5
3:315930000-316250000	.	CDS	84427	84495	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon16
3:315930000-316250000	.	CDS	84559	84658	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon15
3:315930000-316250000	.	CDS	99304	99506	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon14
3:315930000-316250000	.	CDS	100296	100479	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon13
3:315930000-316250000	.	CDS	100611	100673	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon12
3:315930000-316250000	.	CDS	137699	137974	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon11
3:315930000-316250000	.	CDS	150631	150795	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon10
3:315930000-316250000	.	CDS	151035	151179	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon9
3:315930000-316250000	.	CDS	152011	152071	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon8
3:315930000-316250000	.	CDS	152158	152370	.	-	2	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon7
3:315930000-316250000	.	CDS	152470	152561	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon6
3:315930000-316250000	.	CDS	171849	171977	.	-	1	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon5
3:315930000-316250000	.	CDS	199584	199789	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon4
3:315930000-316250000	.	CDS	240355	240510	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon3
3:315930000-316250000	.	CDS	273652	273684	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon2
3:315930000-316250000	.	CDS	313936	314082	.	-	0	ID=CDS:VGSC;Parent=transcript:VGSC_dl;Name=Exon1
###
```
- Example of **ref.mdom.fa** file
A multi alignment FASTA file of VGSC amino-acid (aa) sequences.
The first entry should be house fly (*M. domestica*) VGSC. The second
and third entries are aa sequences of ck and dl transcript variants
of your species, respectively.
This data will be used to convert aa position coordinate of your species to
*M. domestica* universal coordinates.
This file can also be created from **ref.fa** and **ref.bed** using [scripts/make_AA_alignment.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/make_AA_alignment.py)
(MUSCLE and biopython are required).

```
>AAB47604 M.domestica VGSC
MTEDSDSISEEERSLFRPFTRESLLQIEQRIAEHE-KQKELERKRAAEGE----------
QIRYDDEDEDEGPQPDPTLEQGVPIPVRMQGSFPPELASTPLEDIDPFYSNVLTFVVISK
GKDIFRFSASKAMWLLDPFNPIRRVAIYILVHPLFSLFIITTILTNCILMIMPTTPTVES
TEVIFTGIYTFESAVKVMARGFILCPFTYLRDAWNWLDFVVIALAYVTMGIDLGNLAALR
TFRVLRALKTVAIVPGLKTIVGAVIESVKNLRDVIILTMFSLSVFALMGLQIYMGVLTQK
CIKRFPLDGSWGNLTDENWFLHNSNSSNWFTENDGESYPVCGNVSGAGQCGEDYVCLQGF
GPNPNYDYTSFDSFGWAFLSAFRLMTQDFWEDLYQHVLQAAGPWHMLFFIVIIFLGSFYL
VNLILAIVAMSYDELQKKAEEEEAAEEEAIREAEEAAAAKAAKLEERANVAAQAAQDAAD
AAAAALHPEMAKSPT-YSCISYELFVGGEKGNDDNNKEKMSIRSVEVESESVSVIQRQPA
PTTAPATKVRKVSTTSLSLPGSPFNLRRGSRSSHKYTIRNGRGRF-GIPGSDRKPLVLQT
YQDAQQHLPYADDSNAVTPMSEENGAIIVPAYYCNLGSRHSSYTSHQSRISYTSHGDLLG
GMAAMGASTMTKESKLRSRNTRNQSIGAATNGGSSTAGGGYPDANHKEQRDYEMGQDYTD
EAGKIKHHDNPFIEPVQTQTVVDMKDVMVLNDIIEQAAGRHSRASERG------------
---------EDDDEDGPTFKDIALEYILKGIEIFCVWDCCWVWLKFQEWVSFIVFDPFVE
LFITLCIVVNTMFMAMDHHDMNPELEKVLKSGNYFFTATFAIEASMKLMAMSPKYYFQEG
WNIFDFIIVALSLLELGLEGVQGLSVLRSFRLLRVFKLAKSWPTLNLLISIMGRTMGALG
NLTFVLCIIIFIFAVMGMQLFGKNYIDHKDRFKDHELPRWNFTDFMHSFMIVFRVLCGEW
IESMWDCMYVGDVSCIPFFLATVVIGNLVVLNLFLALLLSNFGSSSLSAPTADNDTNKIA
EAFNRIARFKNWVKRNIADCFKLIRNKLTNQISD-QP-----------------------
SEHGDNELELGHDEIMGDGLIKKGMKGETQLEVAIGDGMEFTIHGDMKNNKPKKSKFMNN
TTMIGNSI-NHQDNRLEHELNHRGLSIQDDDTASINSYGSHKNRPFKDESHKGSAETIEG
EEKRDVSKEDLGLDEELDEEAEGDEGQLDGDIIIHAQNDDEIIDDYPADCFPDSYYKKFP
ILAGDEDSPFWQGWGNLRLKTFQLIENKYFETAVITMILMSSLALALEDVHLPDRPVMQD
ILYYMDRIFTVIFFLEMLIKWLALGFKVYFTNAWCWLDFVIVMLSLINLVAVWSGLNDIA
VFRSMRTLRALRPLRAVSRWEGMKVVVNALVQAIPSIFNVLLVCLIFWLIFAIMGVQLFA
GKYFKCKDGNDTVLSHEIIPNRNACKSENYTWENSAMNFDHVGNAYLCLFQVATFKGWIQ
IMNDAIDSREVDKQPIRETNIYMYLYFVFFIIFGSFFTLNLFIGVIIDNFNEQKKKAGGS
LEMFMTEDQKKYYNAMKKMGSKKPLKAIPRPRWRPQAIVFEIVTDKKFDIIIMLFIGLNM
FTMTLDRYDASEAYNNVLDKLNGIFVVIFSGECLLKIFALRYHYFKEPWNLFDVVVVILS
ILGLVLSDIIEKYFVSPTLLRVVRVAKVGRVLRLVKGAKGIRTLLFALAMSLPALFNICL
LLFLVMFIFAIFGMSFFMHVKEKSGINAVYNFKTFGQSMILLFQMSTSAGWDGVLDAIIN
EEDCDPPDNDKGYPGNCGSATVGITFLLSYLVISFLIVINMYIAVILENYSQATEDVQEG
LTDDDYDMYYEIWQQFDPEGTQYIRYDQLSEFLDVLEPPLQIHKPNKYKIISMDMPICRG
DMMYCVDILDALTKDFFARKGNPIEETGEIGEIAARPDTEGYDPVSSTLWRQREEYCAKL
IQNAWRRYKN-------GPPQEGDEGEAAGGEDGAEGGEGEGGSGGGGG-DDGGSATGAT
AAAGATSPSDPDAGEAD-GASVGGPLSPGCVSGGSN---GRQTAVLVESDGFVTKNGHKV
VIHSRSPSITSRTADV
>VGSC_ck
MTEDSDSISEEERSLFRPFTRESLAAIERRIADAEAKQRELEKKR-AEGETGFGRKKKKK
EIRYDDEDEDEGPQPDSTLEQGVPIPVRMQGSFPPELASTPLEDIDSYYANQRTFVVVSK
GKDIFRFSATNALYVLDPFNPIRRVAIYILVHPLFSFFIITTILTNCILMIMPSTPTVES
TEVIFTGIYTFESAVKVMARGFILQPFTYLRDAWNWLDFVVIALAYVTMGIDLGNLAALR
TFRVLRALKTVAIVPGLKTIVGAVIESVKNLRDVIILTMFSLSVFALMGLQIYMGVLTQK
CIREFPMDGSWGNLSDENWERFNNNDSNWYFSETGDT-PLCGNSSGAGQCEEGYICLQGY
GDNPNYGYTSFDTFGWAFLSAFRLMTQDYWENLYQLVLRSAGPWHMLFFIVIIFLGSFYL
VNLILAIVAMSYDELQKKAEEEEAAEEEALREAEEAAAAKAAKLEAQAA-----------
AAAAAANPEIAKSPSDFSCHSYELFVNQEKGNDDNNKEKMSIRSEGLESVSEITRTTAPT
ATAAGTAKARKVSAASLSLPGSPFNLRRGSRGSHQFTIRNGRGRFVGVPGSDRKPLVLST
YLDAQEHLPYADDSNAVTPMSEENGAIIVPVYYANLGSRHSSYTSHQSRISYTSHGDLLG
G--------MTKESRLRNRSARNTNHSIVPPPNMSGPNMSYVDSNHKGQRDFDMSQDCTD
EAGKIKHNDNPFIEPSQTQTVVDMKDVMVLNDIIEQAAGRHSRASDHGVCILSLHRFSFC
CVSVYYFPTEDDDEDGPTFKDKALEFTMRMIDVFCVWDCCWVWLKFQEWVAFIVFDPFVE
LFITLCIVVNTLFMALDHHDMDPDMERALKSGNYFFTATFAIEATMKLIAMSPKYYFQEG
WNIFDFIIVALSLLELGLEGVQGLSVLRSFRLLRVFKLAKSWPTLNLLISIMGRTVGALG
NLTFVLCIIIFIFAVMGMQLFGKNYTDNVDRFPDKDLPRWNFTDFMHSFMIVFRVLCGEW
IESMWDCMLVGDVSCIPFFLATVVIGNLVVLNLFLALLLSNFGSSSLSAPTADNETNKIA
EAFNRISRFSNWIKSNIANALKFVKNKLTSQIASVQPAGEQHNHLSWIWNEGKGVCPCIS
AEHGENELELTPDDILADGLLKKGVKEHNQLEVAIGDGMEFTIHGDLKNKGKKNKQLMNN
SKVIGNSISNHQDNKLEHELNHRGMSLQDDDTASIKSYGSHKNRPFKDESHKGSAETMEG
EEKRDVSKEDLGIDEELDDECDGEEGPLDGELIIHA-DEDEVIEDSPADCCPDNCYKKFP
VLAGDDDAPFWQGWANLRLKTFQLIENKYFETAVITMILLSSLALALEDVHLPHRPILQD
VLYYMDRIFTVIFFLEMLIKWLALGFRVYFTNAWCWLDFIIVMLSLINLTAIWVGAADIP
AFRSMRTLRALRPLRAVSRWEGMRVVVNALVQAIPSIFNVLLVCLIFWLIFAIMGVQLFA
GKYFKCVDKNKTTLSHEIIPDVNACVAENYTWENSPMNFDHVGKAYLCLFQVATFKGWIQ
IMNDAIDSREVGKQPIRETNIYMYLYFVFFIIFGSFFTLNLFIGVIIDNFNEQKKKAGGS
LEMFMTEDQKKYYNAMKKMGSKKPLKAIPRPRWRPQAIVFEIVTNKKFDMIIMLFIGFNM
LTMTLDHYKQTDTFSAVLDYLNMIFICIFSSECLMKIFALRYHYFIEPWNLFDFVVVILS
ILGLVLSDLIEKYFVSPTLLRVVRVAKVGRVLRLVKGAKGIRTLLFALAMSLPALFNICL
LLFLVMFIFAIFGMSFFMHVKDKSGLDDVYNFKTFGQSMILLFQMSTSAGWDGVLDGIIN
EDECLPPDNDKGYPGNCGSATIGITYLLAYLVISFLIVINMYIAVILENYSQATEDVQEG
LTDDDYDMYYEIWQQFDPDGTQYIRYDQLSDFLDVLEPPLQIHKPNKYKIISMDIPICRG
DMMFCVDILDALTKDFFARKGNPIEETAELGEVQARPDEVGYEPVSSTLWRQREEYCARV
IQHAWRKHKERQAGGGGGDDTDADACDNDDGDDGG-GGAGDGGSAGGGGVTSPGVGSGSI
VGGGTTPGSGGGGSQANLGIVVEHNLSPKESPDGNNDPQGRQTAVLVESDGFVTKNGHRV
VIHSRSPSITSRSADV
>VGSC_dl
MTEDSDSISEEERSLFRPFTRESLAAIERRIADAEAKQRELEKKR-AEGETGFGRKKKKK
EIRYDDEDEDEGPQPDSTLEQGVPIPVRMQGSFPPELASTPLEDIDSYYANQRTFVVVSK
GKDIFRFSATNALYVLDPFNPIRRVAIYILVHPLFSFFIITTILTNCILMIMPSTPTVES
TEVIFTGIYTFESAVKVMARGFILQPFTYLRDAWNWLDFVVIALAYVTMGIDLGNLAALR
TFRVLRALKTVAIVPGLKTIVGAVIESVKNLRDVIILTMFSLSVFALMGLQIYMGVLTQK
CIREFPMDGSWGNLSDENWERFNNNDSNWYFSETGDT-PLCGNSSGAGQCEEGYICLQGY
GDNPNYGYTSFDTFGWAFLSAFRLMTQDYWENLYQLVLRSAGPWHMLFFIVIIFLGSFYL
VNLILAIVAMSYDELQKKAEEEEAAEEEALREAEEAAAAKAAKLEAQAA-----------
AAAAAANPEIAKSPSDFSCHSYELFVNQEKGNDDNNKEKMSIRSEGLESVSEITRTTAPT
ATAAGTAKARKVSAASLSLPGSPFNLRRGSRGSHQFTIRNGRGRFVGVPGSDRKPLVLST
YLDAQEHLPYADDSNAVTPMSEENGAIIVPVYYANLGSRHSSYTSHQSRISYTSHGDLLG
G--------MTKESRLRNRSARNTNHSIVPPPNMSGPNMSYVDSNHKGQRDFDMSQDCTD
EAGKIKHNDNPFIEPSQTQTVVDMKDVMVLNDIIEQAAGRHSRASDHGVCILSLHRFSFC
CVSVYYFPTEDDDEDGPTFKDKALEFTMRMIDVFCVWDCCWVWLKFQEWVAFIVFDPFVE
LFITLCIVVNTLFMALDHHDMDPDMERALKSGNYFFTATFAIEATMKLIAMSPKYYFQEG
WNIFDFIIVALSLLELGLEGVQGLSVLRSFRLLRVFKLAKSWPTLNLLISIMGRTMGALG
NLTFVLCIIIFIFAVMGMQLFGKNYIDNVDRFPDKDLPRWNFTDFMHSFMIVFRVLCGEW
IESMWDCMLVGDVSCIPFFLATVVIGNLVVLNLFLALLLSNFGSSSLSAPTADNETNKIA
EAFNRISRFSNWIKSNIANALKFVKNKLTSQIASVQPAGEQHNHLSWIWNEGKGVCPCIS
AEHGENELELTPDDILADGLLKKGVKEHNQLEVAIGDGMEFTIHGDLKNKGKKNKQLMNN
SKVIGNSISNHQDNKLEHELNHRGMSLQDDDTASIKSYGSHKNRPFKDESHKGSAETMEG
EEKRDVSKEDLGIDEELDDECDGEEGPLDGELIIHA-DEDEVIEDSPADCCPDNCYKKFP
VLAGDDDAPFWQGWANLRLKTFQLIENKYFETAVITMILLSSLALALEDVHLPHRPILQD
VLYYMDRIFTVIFFLEMLIKWLALGFRVYFTNAWCWLDFIIVMVSLINFVASLCGAGGIQ
AFKTMRTLRALRPLRAMSRMQGMRVVVNALVQAIPSIFNVLLVCLIFWLIFAIMGVQLFA
GKYFKCVDKNKTTLSHEIIPDVNACVAENYTWENSPMNFDHVGKAYLCLFQVATFKGWIQ
IMNDAIDSREVGKQPIRETNIYMYLYFVFFIIFGSFFTLNLFIGVIIDNFNEQKKKAGGS
LEMFMTEDQKKYYNAMKKMGSKKPLKAIPRPRWRPQAIVFEIVTNKKFDMIIMLFIGFNM
LTMTLDHYKQTDTFSAVLDYLNMIFICIFSSECLMKIFALRYHYFIEPWNLFDFVVVILS
ILGLVLSDLIEKYFVSPTLLRVVRVAKVGRVLRLVKGAKGIRTLLFALAMSLPALFNICL
LLFLVMFIFAIFGMSFFMHVKDKSGLDDVYNFKTFGQSMILLFQMSTSAGWDGVLDGIIN
EDECLPPDNDKGYPGNCGSATIGITYLLAYLVISFLIVINMYIAVILENYSQATEDVQEG
LTDDDYDMYYEIWQQFDPDGTQYIRYDQLSDFLDVLEPPLQIHKPNKYKIISMDIPICRG
DMMFCVDILDALTKDFFARKGNPIEETAELGEVQARPDEVGYEPVSSTLWRQREEYCARV
IQHAWRKHKERQAGGGGGDDTDADACDNDDGDDGG-GGAGDGGSAGGGGVTSPGVGSGSI
VGGGTTPGSGGGGSQANLGIVVEHNLSPKESPDGNNDPQGRQTAVLVESDGFVTKNGHRV
VIHSRSPSITSRSADV
```

### Usage

```
usage: genotyp.py -s species_name -o out_dir_path -l list_file_path
                 (-t num_max_threads -b num_threads_per_bwa
                  -m mode[ngs_dna|ngs_rna] -r ref_dir_path
                  -v variant_caller[freebayes|gatk]
                  )


Genotype VGSC gene.

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
  -v {freebayes,gatk}, --variant_caller {freebayes,gatk}
                        Variant caller to be used. Default is freebayes.
  -n, --no_clean        Do not clean old BAM files after rmdup. Off in
                        default.
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
In default, the **species_dir/** will be searched in the **MoNaS/references/**.
In stead, you can explicitly specify another directory using `-r, --ref_root` option.
If the program could not find **bwadb/**, **hisatdb/**, **ref.fa.fai** or **ref.dict**
within the reference directory, MoNaS will automatically try to create them.

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
#ID             CHROM          POS   REF_ALELE  ALT_ALLELE(s)   GT      QUAL    AA_CHANGE    AA_CHANGE_HOUSEFLY   AD  EXON
Aalb-SP-08      MNAF02001058.1  2207911 T   G  G/G  231859.0   1574F>1574C/1574F>1574C F1534C!!/F1534C!!   0,931   Exon29
Aalb-SP-06      MNAF02001058.1  2207911 T   G  G/G  231859.0   1574F>1574C/1574F>1574C F1534C!!/F1534C!!   0,847   Exon29
Aalb-SP-01      MNAF02001058.1  2207911 T   G  G/G  231859.0   1574F>1574C/1574F>1574C F1534C!!/F1534C!!   1,934   Exon29
Aalb-SP-03      MNAF02001058.1  2207911 T   G  G/G  231859.0   1574F>1574C/1574F>1574C F1534C!!/F1534C!!   0,910   Exon29
Aalb-SP-02      MNAF02001058.1  2207911 T   G  G/G  231859.0   1574F>1574C/1574F>1574C F1534C!!/F1534C!!   1,990   Exon29
Aalb-Toyama-01  MNAF02001058.1  2226519 G   A  G/A  279220.0   wild/2062A>2062T        wild/A2023T         444,472 Exon32
Aalb-Viet-01    MNAF02001058.1  2226519 G   A  G/A  279220.0   wild/2062A>2062T        wild/A2023T         434,372 Exon32
Aalb-SP-08      MNAF02001058.1  2226519 G   A  A/A  279220.0   2062A>2062T/2062A>2062T A2023T/A2023T       0,1048  Exon32
```
Column 1: Sample ID

Column 2: Name of scaffold

Column 3: Nucleotide position of variant

Column 4: Nucleotide(s) in reference

Column 5: Set of alternative alleles

Column 6: Genotype of each sample

Column 7: Quality of variant

Column 8: Annotated AA change **in reference AA position**

Column 9: Annotated AA change **in *M. domestica* AA position**.
          Known *kdr* mutations which are listed in [scripts/kdr_list.json](https://github.com/ItokawaK/MoNaS/blob/master/scripts/kdr_list.json) will be annotated with "!!" character.
          **Note that the list would not be complete!**.

Column 10: Read depth for each allele

Column 11: Corresponding exon

### Aanlyzing Sanger sequence reads
Although it does not seem the best approach, reads from Sanger sequence technology could be analyzed in
MoNaS pipeline by regarding those Sanger reads as NGS reads. [genotype_sanger.py](https://github.com/ItokawaK/MoNaS/blob/master/genotype_sanger.py) is a wrapper script to conduct
chopping input Sanger reads into 150 bp short reads (~ x5 coverage) with fake quality values, writing fastq and sample list files,
and then executing MoNaS for those data.

As any sequecing errors existing in input Sanger reads will be considered as **true variants**, it is important to
trim low-quality regions in advance. Also, ambiguous nucleotide codes (R, Y, S, etc...) are not supported yet (ToDo).

```bash
MoNaS/genotype_sanger.py -s Aalb -t 16 -o out_table.tsv sanger_reads.fa
```


Helper tools
------
1. [scripts/bed2gff.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/bed2gff3.py) creates a gff3 file which is interpretable by bcftools csq from a bed file describing *VGSC* CDSs.
1. [scripts/translate.py](https://github.com/ItokawaK/MoNaS/blob/master/scripts/make_AA_alignment.py) translates genome to VGSC protein using CDS information described in bed file, then conducts
pairwise alignment with VGSC in *M. domestica* which is usable as **ref.mdom.fa**. The script also output mismatched AA
between your reference VGSC and *M. domestica* VGSC as standard error with notation for potential kdr(s) listed in
[scripts/kdr_list.json](https://github.com/ItokawaK/MoNaS/blob/master/scripts/kdr_list.json) if exist in your reference genome.   

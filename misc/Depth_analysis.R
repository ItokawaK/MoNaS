library(tidyverse)

source("create_depth_tbl.R")
source("plot_depth.R")


#Aedes albopictus
## Input parameter
gff_file <- "../MoSo_v0.3/references/gff3/Aalb.gff3"
bam_file_dir <- "../test_dir/BAMs/"
samtools_path <- "/home/itokawa/tools/.INSTALLED/samtools"

sample_names <- c(   #bam files should be [sample_name].bam
  "Aalb-SP-01",
  "Aalb-SP-02",
  "Aalb-SP-03")

depth_tbl_aalb <-
  create_depth_tbl(gff_file, bam_file_dir, sample_names, samtools_path)

plot_depth(depth_tbl_aalb$depth, sample_names[1])

#Aedes aegypti
gff_file <- "../MoSo_v0.3/references/gff3/Aaeg.gff3"
bam_file_dir <- "../20180808_Aaeg/"
samtools_path <- "/home/itokawa/tools/.INSTALLED/samtools"

sample_names <- c(   #bam files should be [sample_name].bam
  "PK69aeg-G-1",
  "PK69aeg-G-2",
  "PK69aeg-G-3")

depth_tbl_aaeg <-
  create_depth_tbl(gff_file, bam_file_dir, sample_names, samtools_path)

plot_depth(depth_tbl_aaeg$depth, sample_names[1])


#Culex quinquefasciatus

gff_file <- "../MoSo_v0.3/references/gff3/Cpip.gff3"
bam_file_dir <- "../201712_Cpip/BAMs/"
samtools_path <- "/home/itokawa/tools/.INSTALLED/samtools"

sample_names <- c(   #bam files should be [sample_name].bam
  "JPP-01",
  "JPP-02",
  "JNA-01",
  "JNA-02")

depth_tbl_cpip <-
  create_depth_tbl(gff_file, bam_file_dir, sample_names, samtools_path)

plot_depth(depth_tbl_cpip$depth, sample_names[1])
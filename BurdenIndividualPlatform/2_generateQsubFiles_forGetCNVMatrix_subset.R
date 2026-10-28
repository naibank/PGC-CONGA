# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if(!dir.exists("qsubCNVMatrix"))
  dir.create("qsubCNVMatrix")

if(!dir.exists("cnvMatrix"))
  dir.create("cnvMatrix")

args = commandArgs(trailingOnly=TRUE)
gs_path <- args[1]
sample_path <- args[2]
cnv_path <- args[3]
genesets_per_file <- as.numeric(args[4])
queue <- args[5]
# message(gs_path)

disorder <- gsub("sampleinfo|\\.|tsv", "", basename(sample_path))

gs_variable <- ls()

load(gs_path)
gs_variable <- setdiff(ls(), c(gs_variable, "gs_variable"))
gsData.s <- get(gs_variable)

if(basename(gs_path) == "human_gsData_goincludingIEA100-1000_pathways15-500_20211125.RData"){
  gsData.s <- gsData.s$gs2gene
  gsData.s <- gsData.s[sapply(gsData.s, length) >= 50 & sapply(gsData.s, length) <= 500]

  names(gsData.s) <- gsub(":", "_", names(gsData.s))
  names(gsData.s) <- gsub("-", "_", names(gsData.s))
  names(gsData.s) <- gsub("\\.", "_", names(gsData.s))
}
names(gsData.s) <- gsub("\\ |-", ".", names(gsData.s))

nsub <- genesets_per_file
files <- seq(1, length(gsData.s), by = nsub)
sh <- c()
qsub <- readLines(sprintf("/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/main_scripts/GenesetBased_fast/run.parallel.template.%s.sh", queue))
for(f in files){
  fend <- ifelse(f + nsub - 1 < length(gsData.s), f+nsub - 1, length(gsData.s))
  qsub[11] <- sprintf("Rscript /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/main_scripts/GenesetBased_fast/getCNVMatrix_subset.R --gs_path %s --sample_path %s --cnv_path %s --from %s --end %s",
                      gs_path, sample_path, cnv_path, f, fend)
  qsub[3] <- sprintf("#SBATCH --job-name %s.%s", disorder, f)
  sh <- c(sh, sprintf("sbatch qsub.%s.sh", f))
  writeLines(qsub, sprintf("qsubCNVMatrix/qsub.%s.sh", f))
}

writeLines(sh, "qsubCNVMatrix/run.sh")

system("cd qsubCNVMatrix; sh run.sh")

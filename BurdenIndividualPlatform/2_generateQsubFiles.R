# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if(!dir.exists("qsub"))
  dir.create("qsub")

if(!dir.exists("results"))
  dir.create("results")

args = commandArgs(trailingOnly=TRUE)
gs_path <- args[1]
disorder <- args[2]
queue <- args[3]

files <- list.files("cnvMatrix/", full.names = T)
# queue <- "tcag"
sh <- c()
qsub <- readLines(sprintf("/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/main_scripts/GenesetBased_fast/run.parallel.template.%s.sh", queue))
for(f in files){
  qsub[11] <- sprintf("Rscript /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/main_scripts/GenesetBased_fast/burden.R --disorder %s --file ../%s --gs_path %s",
                      disorder, f, gs_path)
  filename <- gsub("cnvMatrix\\/\\/|.tsv", "", basename(f))
  qsub[3] <- sprintf("#SBATCH --job-name %s", filename)
  sh <- c(sh, sprintf("sbatch qsub.%s.sh", filename))
  writeLines(qsub, sprintf("qsub/qsub.%s.sh", filename))
}

writeLines(sh, "qsub/run.sh")

system("cd qsub; sh run.sh")

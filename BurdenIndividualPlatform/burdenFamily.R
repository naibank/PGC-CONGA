#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[1]))), "fast_helpers.R"))

# devtools::install_github("naibank/GSBurden")
library(data.table)
library(survival)
library(foreach)
library(doParallel)

# One thread per worker avoids nested data.table/BLAS parallelism.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
data.table::setDTthreads(1L)

CNVBurdenTest <- function (cnv.matrix, geneset, label, covariates, correctGlobalBurden = T, 
                           standardizeCoefficient = T) 
{
  distinct.prefixes <- c("DEL", "DUP", "DELDUP")
  if (!all(cnv.matrix[[label]][!is.na(cnv.matrix[[label]])] %in% c(0, 1)) ||
      length(unique(na.omit(cnv.matrix[[label]]))) != 2L)
    stop("This script requires a binary outcome coded 0/1 with both groups present")
  required <- unique(c(label, covariates, "FID", "gene_count_DEL", "gene_count_DUP",
                       paste0(names(geneset), "_DEL"), paste0(names(geneset), "_DUP")))
  if (length(setdiff(required, names(cnv.matrix))))
    stop("Missing matrix columns: ", paste(setdiff(required, names(cnv.matrix)), collapse=", "))

  cnv.matrix <- cnv.matrix[, required, drop=FALSE]
  test.out <- data.frame()
  for (cnvtype in distinct.prefixes) {
    out.message <- sprintf("Testing for %s CNVs", cnvtype)
    message(out.message)
    
    test.tmp <- data.frame()
    test.tmp <- foreach (this.gs = names(geneset), .combine=rbind) %dopar% {
      library(survival)
      model.data <- cnv.matrix[, unique(c(label, covariates, "FID")), drop=FALSE]
      if(cnvtype == "DELDUP"){
        # feature <- sprintf("%s_%s", this.gs, cnvtype)
        # global <- sprintf("gene_count_%s", cnvtype)
        this.covariates <- covariates
        
        model.data$feature <- cnv.matrix[, sprintf("%s_%s", this.gs, "DEL")] + cnv.matrix[, sprintf("%s_%s", this.gs, "DUP")]
        model.data$global <- cnv.matrix[, sprintf("gene_count_%s", "DEL")] + cnv.matrix[, sprintf("gene_count_%s", "DUP")] -  model.data$feature
        
      }else{
        feature <- sprintf("%s_%s", this.gs, cnvtype)
        global <- sprintf("gene_count_%s", cnvtype)
        this.covariates <- covariates
        
        model.data$feature <- cnv.matrix[, feature]
        model.data$global <- cnv.matrix[, global] -  cnv.matrix[, feature]
      }

      if (!correctGlobalBurden) model.data$global <- NULL
      model.data <- model.data[complete.cases(model.data), , drop=FALSE]
      feature.sd <- sd(model.data$feature)
      if (!nrow(model.data) || !is.finite(feature.sd) || feature.sd == 0) {
        return(data.frame(geneset=this.gs, type=cnvtype, coefficient=NA_real_,
                          coeff.upper=NA_real_, coeff.lower=NA_real_, pvalue=NA_real_,
                          stderr=NA_real_, testval=NA_real_, waldp=NA_real_))
      }
      if (standardizeCoefficient)
        model.data$feature <- as.numeric(scale(model.data$feature))
      
      ref.term <- reformulate(c(covariates, "strata(FID)", if (correctGlobalBurden) "global"), response=label)
      ref.model <- clogit(ref.term, model.data)
      add.model <- clogit(update(ref.term, . ~ . + feature), model.data)
      ano <- anova(ref.model, add.model, test = "Chisq")
      
      names(ano)[length(names(ano))] <- "pvalue"
      pvalue <- ano$pvalue[2]
      coefficient <- add.model$coefficients["feature"]
      sm <- summary(add.model)
      if ("feature" %in% rownames(sm$coefficients)) {
        waldp <- sm$coefficients["feature", 5]
        stderr <- sm$coefficients["feature", 3]
        testval <- sm$coefficients["feature", 4]
      }
      else {
        waldp <- stderr <- testval <- NA
      }
      
      conf <- tryCatch({
        confint(add.model)
      }, error = function(e) {
        return(NA)
      })
      if (is.null(nrow(conf))) {
        coeff.l <- NA_real_
        coeff.u <- NA_real_
      }
      else {
        coeff.l <- conf["feature", 1]
        coeff.u <- conf["feature", 2]
      }
      temp.out <- data.frame(geneset = this.gs, type = cnvtype, 
                             coefficient = coefficient, coeff.upper = coeff.u, 
                             coeff.lower = coeff.l, pvalue = pvalue, stderr, 
                             testval, waldp)
      temp.out
    }
    test.out <- rbind(test.out, test.tmp)
  }
  
  test.out <- test.out[order(test.out$pvalue), ]

  return(test.out)
}

#### read arguments from command line
params <- parse_args(commandArgs(trailingOnly=TRUE), c("gs_path", "disorder", "file"))
# --set all --from 1 --end 12 --gs_path /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenAnalysis/geneset_celltype/gsCelltype.RData --disorder asd_fam
### params set, from, end
gs_path <- params$gs_path
disorder <- params$disorder
infile <- params$file
filename <- sub("\\.tsv$", "", basename(infile))
filename <- gsub(paste0(disorder, "\\.matrix\\."), "", filename)
#############################
#############################
gsData.s <- load_genesets(gs_path)

cnv.matrix <- fread(infile, data.table = F)
dt.colsum <- colSums(cnv.matrix[, setdiff(grep("DEL|DUP", names(cnv.matrix)),
                                          grep("cnv_count|cnv_size|gene_count", names(cnv.matrix))), drop=FALSE] != 0, na.rm=TRUE)

dt.colsum <- data.frame("geneset" = names(cnv.matrix)[setdiff(grep("DEL|DUP", names(cnv.matrix)),
                                                              grep("cnv_count|cnv_size|gene_count", names(cnv.matrix)))],
                        "samplecount" = dt.colsum)

{
  dir.create("../colSums", showWarnings=FALSE, recursive=TRUE)
  write.table(dt.colsum, sprintf("../colSums/%s.%s.burden.tsv", disorder, filename), sep="\t", row.names=F, quote=F, col.names=T)
}

covariates <- c("sex", paste0("PC", c(1:10)))
if(length(unique(cnv.matrix$PLATFORM)) == 1){
  covariates <- c("sex", paste0("PC", c(1:10)))
}
#### perform global burden and gene set burden test ####
from <- strsplit(filename, "\\.")[[1]][2]
to <- strsplit(filename, "\\.")[[1]][3]
gsData.s <- subset_genesets(gsData.s, from, to)

run_tests <- function() {
  ncores <- suppressWarnings(as.integer(Sys.getenv("BURDEN_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))
  if (is.na(ncores) || ncores < 1L) stop("BURDEN_CORES must be a positive integer")
  ncores <- min(ncores, length(gsData.s))
  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    on.exit({parallel::stopCluster(cl); foreach::registerDoSEQ()}, add=TRUE)
    doParallel::registerDoParallel(cl)
  } else foreach::registerDoSEQ()
if(length(gsData.s) > 0){
  for(platform in unique(cnv.matrix$PLATFORM)){
    dt <- CNVBurdenTest(cnv.matrix[which(cnv.matrix$PLATFORM == platform), ], gsData.s, "status", covariates)
    
    if(!dir.exists(sprintf("../results/%s/%s", disorder, 
                           platform))){
      dir.create(sprintf("../results/%s/%s", disorder, 
                         platform), recursive = T)
    }
    
    #### write output to files ####
    
    write.table(dt, sprintf("../results/%s/%s/%s.%s.burden.tsv", disorder, platform, disorder, filename), sep="\t", row.names=F, quote=F, col.names=T)
  }
}

}
run_tests()

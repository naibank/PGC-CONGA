#########################################################
#### Install GSBurden ####
### only run a line below when doing this in Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
source(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[1]))), "fast_helpers.R"))

# devtools::install_github("naibank/GSBurden")
# Rscript 1_getCNVMatrix_ASD.R <geneset_path> <sample_info_path> <cnv_path>
library(data.table)
library(GenomicRanges)
library(GSBurden)

path <- "/hpf/largeprojects/"

if(!dir.exists(path)){
  path <- "/Volumes/"
}

params <- parse_args(commandArgs(trailingOnly=TRUE), c("gs_path", "sample_path", "cnv_path", "from", "end"))
# --set all --from 1 --end 12 --gs_path /hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenAnalysis/geneset_celltype/gsCelltype.RData --disorder asd
### params set, from, end
gs_path <- params$gs_path
sample_path <- params$sample_path
cnv_path <- params$cnv_path
from <- params$from
end <- params$end


### for testing
# gs_path <- "/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/geneset_celltype/gsCelltype.RData"
# sample_path <- "/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/BurdenIndividualPlatform/genebased_mega_analysis/BD/sampleinfo.bd.tsv"
# cnv_path <- "/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/Data/CNVPlinkFormat/ALL/BD/noFam/all_platforms_10kb_10probes_noFams_BD_ALL.cnv"
# gene_level_sumstat <- "/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2024/MetaAnalysis/genebased_mega_analysis/BD/BD_final_summ_stat.tsv"

disorder <- gsub("sampleinfo|\\.|tsv", "", basename(sample_path))

gene.in <- read.delim(paste0(path, "/tcagstor/users/worrawat/PGC/2026/FromOmar/FilteredExons/hg38_refGene_20200708_noExonWith20CNVFiltered.txt"), stringsAsFactors = F)
genes <- GeneAnnotation(gene.in$enzid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)
genes.g <- GRanges(genes$chr, IRanges(genes$start, genes$end), "*")

gsData.s <- load_genesets(gs_path)
gsData.s <- subset_genesets(gsData.s, from, end)

sampleinfo <- data.table::fread(sample_path, data.table = F)
if (ncol(sampleinfo) != 18L) stop("Expected 18 sample-info columns")
if (anyNA(sampleinfo[[2]]) || anyDuplicated(sampleinfo[[2]])) stop("Sample IDs must be nonmissing and unique")
names(sampleinfo) <- c("FID", "sample", "status", "sex", "ancestry", "PLATFORM", "dataset", "phenotype", paste0("PC", 1:10))

cnvs.in <- fread(cnv_path, data.table = F)
if (ncol(cnvs.in) != 8L) stop("Expected eight CNV columns")
names(cnvs.in) <- c("FID", "sample", "chr", "start", "end", "type", "score", "site")
cnvs.in <- cnvs.in[which(cnvs.in$type %in% c(1, 3)), ]
cnvs.in$chr <- sub("^23$", "X", sub("^chr", "", cnvs.in$chr))
cnvs.in$chr <- sub("^24$", "Y", cnvs.in$chr)
cnvs.in$chr <- paste0("chr", cnvs.in$chr)
cnvs.in$type <- ifelse(cnvs.in$type == 1, "DEL", "DUP")

cnvs.in <- cnvs.in[cnvs.in$sample %in% sampleinfo$sample, ]
cnvs.in <- cnvs.in[cnvs.in$chr %in% paste0("chr", c(1:22, "X")),]
cnvs.in$cnvid <- paste(cnvs.in$chr, cnvs.in$start, cnvs.in$end, cnvs.in$type, sep="#")

cnvs.in.g <- GRanges(cnvs.in$chr, IRanges(cnvs.in$start, cnvs.in$end), "*")
olap <- data.frame(findOverlaps(cnvs.in.g, genes.g))
olap$gene <- genes$enzid[olap$subjectHits]
nrxn1_del <- unique(olap$queryHits[genes$gsymbol[olap$subjectHits] == "NRXN1" &
                              cnvs.in$type[olap$queryHits] == "DEL"])
nrxn1_del <- unique(cnvs.in$sample[nrxn1_del])

# Reuse the overlap calculation for both types and every geneset.
# Count unique sample/gene pairs, matching GSBurden::getCNVGSMatrix.
build_matrix <- function(type) {
  ids <- unique(cnvs.in$sample)
  out <- data.frame(sample=ids, stringsAsFactors=FALSE)
  rows <- which(cnvs.in$type == type)
  sample.index <- match(cnvs.in$sample[rows], ids)
  out[[paste0("cnv_count_", type)]] <- tabulate(sample.index, nbins=length(ids))
  sizes <- cnvs.in$end[rows] - cnvs.in$start[rows] + 1
  totals <- if (length(sample.index)) rowsum(matrix(sizes, ncol=1), sample.index, reorder=FALSE) else NULL
  size.out <- numeric(length(ids))
  if (length(sample.index)) size.out[as.integer(rownames(totals))] <- totals[, 1]
  out[[paste0("cnv_size_", type)]] <- size.out
  hits <- which(cnvs.in$type[olap$queryHits] == type)
  pairs <- unique(data.frame(sample=match(cnvs.in$sample[olap$queryHits[hits]], ids),
                             gene=olap$gene[hits]))
  out[[paste0("gene_count_", type)]] <- tabulate(pairs$sample, nbins=length(ids))
  for (gs in names(gsData.s))
    out[[paste0(gs, "_", type)]] <- tabulate(pairs$sample[pairs$gene %in% gsData.s[[gs]]], nbins=length(ids))
  out
}
del.matrix <- build_matrix("DEL")
dup.matrix <- build_matrix("DUP")
cnv.matrix <- cbind(del.matrix, dup.matrix[, setdiff(names(dup.matrix), "sample"), drop=FALSE])
rm(del.matrix, dup.matrix)

platform.case <- table(sampleinfo$PLATFORM[sampleinfo$status == 1])
platform.case <- names(which(platform.case >= 50))
platform.ctrl <- table(sampleinfo$PLATFORM[sampleinfo$status == 0])
platform.ctrl <- names(which(platform.ctrl >= 50))
platform.both <- intersect(platform.case, platform.ctrl)

sampleinfo <- sampleinfo[sampleinfo$PLATFORM %in% platform.both,]

cnv.matrix <- cnv.matrix[cnv.matrix$sample %in% sampleinfo$sample, ]
cnv.matrix <- merge(sampleinfo, cnv.matrix, by = "sample",
                    all.x = T)
# Zero-fill burden columns only: missing phenotype/covariates remain missing.
for (column in setdiff(names(cnv.matrix), names(sampleinfo)))
  cnv.matrix[[column]][is.na(cnv.matrix[[column]])] <- 0

dir.create("../cnvMatrix", recursive=TRUE, showWarnings=FALSE)
write.table(cnv.matrix, sprintf("../cnvMatrix/%s.matrix.all.%s.%s.tsv", disorder, from, end), sep="\t", row.names=F, quote=F, col.names=T)

### get cnv.matrix excluding recurrent CNVs + genome-wide significant genes
recurrent_loci <- fread("../../../../../2026/FromOmar/GWSLoci/latest_gws_cnvs.txt", data.table = F)
recurrent_loci <- recurrent_loci[recurrent_loci$NAME != "NRXN1", ]
recurrent_loci.g <- GRanges(paste0("chr", recurrent_loci$CHR), IRanges(as.numeric(recurrent_loci$START), as.numeric(recurrent_loci$END)), "*")

recurrent_olap <- data.frame(findOverlaps(cnvs.in.g, recurrent_loci.g))
recurrent_olap$coverage <- width(pintersect(cnvs.in.g[recurrent_olap$queryHits],
                                            recurrent_loci.g[recurrent_olap$subjectHits]))/width(recurrent_loci.g[recurrent_olap$subjectHits])
recurrent_olap$sample <- cnvs.in$sample[recurrent_olap$queryHits]
recurrent_olap$type <- cnvs.in$type[recurrent_olap$queryHits]
# Empty overlaps must still produce a valid empty result.
recurrent_coverage <- as.data.table(recurrent_olap)[,
  list(coverage=sum(coverage)), by=c("sample", "type", "subjectHits")]
recurrent_coverage <- as.data.frame(recurrent_coverage)
recurrent_coverage <- recurrent_coverage[which(recurrent_coverage$coverage > 0.75), ]
recurrent_coverage$locus <- recurrent_loci$locus_name[recurrent_coverage$subjectHits]
recurrent_coverage <- merge(recurrent_coverage, cnv.matrix[, c("sample", "status")], by = "sample", all.x = T)

excluding.samples <- unique(c(nrxn1_del, recurrent_coverage$sample))

if(as.integer(from) == 1){
  writeLines(excluding.samples, "../samples.excluded.with.known.significant.variants.txt")
}

cnv.matrix <- cnv.matrix[which(!cnv.matrix$sample %in% excluding.samples), ]
write.table(cnv.matrix, sprintf("../cnvMatrix/%s.matrix.smaller.%s.%s.tsv", disorder, from, end), sep="\t", row.names=F, quote=F, col.names=T)


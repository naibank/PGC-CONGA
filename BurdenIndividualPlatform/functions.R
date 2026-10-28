combineResults <- function(disorder, set = "all", manifest, samplewithbigcnvs, cnv, genes, gsData.s, gs, dropMismatch=F){
   
    platform.case <- table(manifest$PLATFORM[manifest$status == 1])
    platform.case <- names(which(platform.case >= 100))
    platform.ctrl <- table(manifest$PLATFORM[manifest$status == 0])
    platform.ctrl <- names(which(platform.ctrl >= 100))
    platform.both <- intersect(platform.case, platform.ctrl)
    
    if(set == "smaller")
      manifest <- manifest[!manifest$sample %in% samplewithbigcnvs, ]
    
    manifest <- manifest[manifest$PLATFORM %in% platform.both,]
    
    names(cnv) <- c("FID", "sample", "chr", "start", "end", "type", "score", "site")
    cnv <- cnv[which(cnv$type %in% c(1,3)), ]
    cnv$chr <- gsub("23", "X", cnv$chr)
    cnv$chr <- gsub("24", "Y", cnv$chr)
    cnv$chr <- paste0("chr", cnv$chr)
    cnv$type <- ifelse(cnv$type == 1, "DEL", "DUP")
    cnv <- cnv[cnv$sample %in% manifest$sample, ]
    cnv <- merge(cnv, manifest[, c("sample", "status")], by="sample", all.x = T)
    names(genes) <- c("chr", "start", "end", "isoform", "genesymbol", "enzid")
    genes.g <- GRanges(genes$chr, IRanges(genes$start, genes$end), '*')
    cnv.g <- GRanges(cnv$chr, IRanges(cnv$start, cnv$end), "*")
    
    olap <- data.frame(findOverlaps(cnv.g, genes.g))
    olap$sample <- cnv$sample[olap$queryHits]
    olap$status <- cnv$status[olap$queryHits]
    olap$gene <- genes$genesymbol[olap$subjectHits]
    olap$type <- cnv$type[olap$queryHits]
    olap <- unique(olap[, c("sample", "status", "gene", "type")])
    olap <- as.data.frame.matrix(table(paste(olap$gene, olap$type, sep="#"), olap$status))
    olap$gene <- sapply(sapply(rownames(olap), strsplit, "#"), "[", 1) 
    olap$type <- sapply(sapply(rownames(olap), strsplit, "#"), "[", 2)
    olap <- olap[, c(3, 4, 1,2)]
    names(olap) <- c("gene", "type", "control", "case")
    
    # gs <- gsData.s$gs2gene
    # names(gs) <- gsub(":|\\.|-", "_", names(gs))
    # 
    # gsData.s <- gsData.s$gs2name
    # names(gsData.s) <- gsub(":|\\.|-", "_", names(gsData.s))
    # gsData.s <- data.frame("geneset" = names(gsData.s), "fullname" = gsData.s)
    ncases <- sum(manifest$status == 1)
    nctrls <- sum(manifest$status == 0)
    olap$caseFreq <- olap$case/ncases
    olap$ctrlFreq <- olap$control/nctrls
    
    write.table(olap, sprintf("%s.gene.freq.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    olap$exp <- as.numeric(olap$caseFreq > olap$ctrlFreq)
    olap$exp[which(olap$caseFreq < olap$ctrlFreq)] <- -1
    
    write.table(olap[olap$type == "DEL", c("gene", "exp")], sprintf("%s.del.exp.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    write.table(olap[olap$type == "DUP", c("gene", "exp")], sprintf("%s.dup.exp.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    
    golap <- merge(olap, unique(genes[, c("genesymbol", "enzid")]), by.x = "gene", by.y = "genesymbol", all.x = T)
    remove.gs <- readLines("/hpf/largeprojects/tcagstor/users/worrawat/PGC/freeze_2022/BurdenAnalysis/bak_2024/pre_2024_backups/AnalysisVariations/geneset_pipeline/GO_to_drop.txt") # redundant GO with 100% overlapped genes
    remove.gs <- gsub(":", "_", remove.gs)
    remove.gs <- gsub("-", "_", remove.gs)
    
    # for(set in c("all", "smaller")){
    files <- list.files("results/", full.names = T)
    tmp.files <- files[grep(sprintf("\\.%s\\.",set), files)]  
    
    dt <- data.frame()
    for(f in tmp.files){
      tmp <- read.delim(f, stringsAsFactors = F)
      dt <- rbind(dt, tmp)
    }
    
    dt <- merge(dt, gsData.s, by = "geneset", all.x = T)
    dt$phenotype <- sign(dt$coefficient)
    dt <- dt[which(!is.na(dt$coefficient) & !is.na(dt$pvalue)), ]
    dt <- dt[which(!dt$geneset %in% remove.gs), ]
    
    dt$qvalue <- NA
    dt$qvalue[dt$type == "DEL"] <- p.adjust(dt$pvalue[dt$type == "DEL"], method = "BH")
    dt$qvalue[dt$type == "DUP"] <- p.adjust(dt$pvalue[dt$type == "DUP"], method = "BH")
    dt$qvalue[dt$type == "All"] <- p.adjust(dt$pvalue[dt$type == "All"], method = "BH")
    
    dt$case <- 0
    dt$control <- 0
    dt$OR <- NA
    dt$genes <- ""
    
    message(nrow(dt))
    # dt <- dt[!dt$geneset %in% remove.gs, ]
    # dt <- dt[grep("^GO", dt$geneset), ]
    
    for(i in 1:nrow(dt)){
      type <- strsplit(dt$type[i], "_")[[1]][1]
      if(type == "All"){
        type <- c("DEL", "DUP")
      }
      
      count <- golap[golap$enzid %in% gs[[dt$geneset[i]]] &
                       golap$type %in% type, ]
      # denom.count <- golap[!golap$enzid %in% gs[[dt$geneset[i]]] &
      #                        golap$type == dt$type[i], ]
      
      dt$case[i] <- sum(count$case)
      dt$control[i] <- sum(count$control)
      
      # denom.case <- sum(denom.count$case)
      # denom.control <- sum(denom.count$control)
        
      dt$OR[i] <- (dt$case[i]/ncases)/(dt$control[i]/nctrls)
      #message(sprintf("OR=%s, coeff=%s, case=%s, ctrl=%s, ncases=%s, nctrl=%s, type=%s", dt$OR[i], dt$coefficient[i], dt$case[i], dt$control[i], ncases, nctrls, dt$type[i]))
      
      if(dt$OR[i] > 1){
        count <- count[count$exp == 1, ]
      }else{
       	count <- count[count$exp == -1, ]
      }
      
      
      dt$genes[i] <- paste(count$gene, collapse = ",")
    }
    
    if(dropMismatch)
    dt <- dt[(dt$OR < 1 & dt$phenotype == -1) |
               (dt$OR > 1 & dt$phenotype == 1), ]
    
    write.table(dt, sprintf("%s.gs.result.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    
    dt <- dt[dt$qvalue < 0.05,]
    
    del <- dt[dt$type == "DEL", ]
    dup <- dt[dt$type == "DUP", ]
    
    write.table(dt, sprintf("%s.gs.result.q0.05.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    
    write.table(del[, c("fullname", "geneset", "pvalue", "qvalue", "phenotype", "genes")], 
                sprintf("%s.del.gs.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
    
    write.table(dup[, c("fullname", "geneset", "pvalue", "qvalue", "phenotype", "genes")], 
                sprintf("%s.dup.gs.%s.tsv", tolower(disorder), set), sep="\t", row.names=F, quote=F, col.names=T)
}

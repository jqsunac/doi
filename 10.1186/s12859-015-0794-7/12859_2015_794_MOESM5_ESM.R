## R code for obtaining results of Blekhman’s count data. After execution of
## this R-code, full results of real data analysis can be obtained. 


library(TCC)            # v1.7.15
library(samr)           # v2.0
library(PoissonSeq)     # v1.1.2
library(EBSeq)          # v1.6.0
library(limma)          # v3.22.2
library(VennDiagram)    # 1.6.9
library(gclus)          # 1.3.1


##---------------------------------------------------------------------------##
## paramter                                                                  ##
##---------------------------------------------------------------------------##
#in_f <- "http://genome.cshlp.org/content/suppl/2009/12/16/gr.099226.109.DC1/suppTable1.xls"
in_f <- "suppTable1.xls"  # input file
param_FDR <- 0.05         # FDR shreshold
methods <- c("EEE-E", "DDD-D", "SSS-S", "E-E(edgeR)", "edgeR_robust", "D-D(DESeq)", 
             "S-S(DESeq2)", "voom", "SAMseq", "PoissonSeq", "baySeq", "EBSeq")

couple <- matrix(1:6, nrow=3, byrow=T)

##---------------------------------------------------------------------------##
## data                                                                      ##
##---------------------------------------------------------------------------##
dataOri <- read.table(in_f, header=TRUE)
data36 <- dataOri[,21:56]
rownames(data36) <- dataOri[, 1]

Ngene <- nrow(dataOri)
data.ori <- cbind(
  dataOri$R1L4.HSF1 + dataOri$R4L2.HSF1, dataOri$R2L7.HSF2 + dataOri$R3L2.HSF2,
  dataOri$R8L1.HSF3 + dataOri$R8L2.HSF3, dataOri$R1L1.HSM1 + dataOri$R5L2.HSM1,
  dataOri$R2L3.HSM2 + dataOri$R4L8.HSM2, dataOri$R3L6.HSM3 + dataOri$R4L1.HSM3,
  dataOri$R1L2.PTF1 + dataOri$R4L4.PTF1, dataOri$R2L4.PTF2 + dataOri$R6L6.PTF2,
  dataOri$R3L7.PTF3 + dataOri$R5L3.PTF3, dataOri$R1L6.PTM1 + dataOri$R3L3.PTM1,
  dataOri$R2L8.PTM2 + dataOri$R4L6.PTM2, dataOri$R6L2.PTM3 + dataOri$R6L4.PTM3,
  dataOri$R1L7.RMF1 + dataOri$R5L1.RMF1, dataOri$R2L2.RMF2 + dataOri$R5L8.RMF2,
  dataOri$R3L4.RMF3 + dataOri$R4L7.RMF3, dataOri$R1L3.RMM1 + dataOri$R3L8.RMM1,
  dataOri$R2L6.RMM2 + dataOri$R5L4.RMM2, dataOri$R3L1.RMM3 + dataOri$R4L3.RMM3)
rownames(data.ori) <- dataOri[, 1]
colnames(data.ori) <- c(
  "HS_rep1", "HS_rep2", "HS_rep3", "HS_rep4", "HS_rep5", "HS_rep6", 
  "PT_rep1", "PT_rep2", "PT_rep3", "PT_rep4", "PT_rep5", "PT_rep6", 
  "RM_rep1", "RM_rep2", "RM_rep3", "RM_rep4", "RM_rep5", "RM_rep6")
HS <- data.ori[, 1:6]
PT <- data.ori[, 7:12]
RM <- data.ori[, 13:18]

data18 <- data.ori
colnames(data18) <- c(
  "HSF1", "HSF2", "HSF3", "HSM1", "HSM2", "HSM3", 
  "PTF1", "PTF2", "PTF3", "PTM1", "PTM2", "PTM3", 
  "RMF1", "RMF2", "RMF3", "RMM1", "RMM2", "RMM3")

##---------------------------------------------------------------------------##
## Sample clustering of raw count data                                       ##
##---------------------------------------------------------------------------##
png(file = "Additional6a.png", width = 800, height = 350)
hc <- clusterSample(data36, dist.method = "spearman", 
             hclust.method = "average", unique.pattern = TRUE)
par(mar=c(0, 4, 0, 0))
plot(hc, sub="", xlab ="", ylab="Height", main="", cex=1.3, cex.lab=1.2)
dev.off()

png(file = "Additional6b.png", width = 800, height = 350)
hc <- clusterSample(data18, dist.method = "spearman", 
             hclust.method = "average", unique.pattern = TRUE)
par(mar=c(0, 4, 0, 0))
plot(hc, sub="", xlab ="", ylab="Height", main="", cex=1.3, cex.lab=1.2)
dev.off()

png(file = "Additional6c.png", width = 800, height = 350)
hc <- clusterSample(data.ori, dist.method = "spearman", 
             hclust.method = "average", unique.pattern = TRUE)
par(mar=c(0, 4, 0, 0))
plot(hc, sub="", xlab ="", ylab="Height", main="", cex=1.3, cex.lab=1.2)
dev.off()


##---------------------------------------------------------------------------##
## Main process                                                              ##
##---------------------------------------------------------------------------##
list <- list()
for(m in 1:length(methods)) {
  if(length(grep(methods[m], names(list))) == 0) {
    list[[methods[m]]] <- matrix(0, nrow = Ngene, ncol = 4)
  } 
}

for(p in 1:4) {
  if (p == 4) {
    Nrep <- 6
    data <- data.ori
  }else{
    Nrep <- 2
    data <- cbind(HS[, couple[p, ]], PT[, couple[p, ]], RM[, couple[p, ]])
  }
  group <- gl(3, Nrep)
  design <- model.matrix(~ group)
  tcc <- new("TCC", data, group)

  ##  EEE-E
  eeee <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
  DE <- estimateDE(eeee, test.method = "edger", design = design, coef = 2:3)
  res <- getResult(DE)
  FDR <- res$q.value
  rank <- res$rank
  res.eeee <- as.data.frame(cbind(FDR, rank))
  rownames(res.eeee) <- rownames(data)
  degs.eeee <- rownames(data[FDR < param_FDR, ])
  list[["EEE-E"]][ , p] <- rownames(res.eeee[order(res.eeee$rank), ])

  ##  DDD-D
  dddd <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", iteration = 3)
  DE <- estimateDE(dddd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
  res <- getResult(DE)
  FDR <- res$q.value
  rank <- res$rank
  res.dddd <- as.data.frame(cbind(FDR, rank))
  rownames(res.dddd) <- rownames(data)
  degs.dddd <- rownames(data[FDR < param_FDR, ])
  list[["DDD-D"]][ , p] <- rownames(res.dddd[order(res.dddd$rank), ])

  ##  SSS-S
  ssss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2", iteration = 3)
  DE <- estimateDE(ssss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
  res <- getResult(DE)
  FDR <- res$q.value
  rank <- res$rank
  res.ssss <- as.data.frame(cbind(FDR, rank))
  rownames(res.ssss) <- rownames(data)
  degs.ssss <- rownames(data[FDR < param_FDR, ])
  list[["SSS-S"]][ , p] <- rownames(res.ssss[order(res.ssss$rank), ])

  ##  E-E(edgeR)
  #ee <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
  #DE <- estimateDE(ee, test.method = "edger", design = design, coef = 2:3)
  ee <- DGEList(counts = data, group = group)
  ee <- edgeR::calcNormFactors(ee, method = "TMM")
  ee <- edgeR::estimateGLMCommonDisp(ee, design, method = "CoxReid")
  ee <- edgeR::estimateGLMTrendedDisp(ee, design, method = "auto")
  ee <- edgeR::estimateGLMTagwiseDisp(ee, design)
  ee <- edgeR::glmFit(ee, design)
  ee <- edgeR::glmLRT(ee, coef = 2:3, test = "chisq")
  res <- ee$table
  pval <- res$PValue
  pval[is.na(pval)] <- 1
  FDR <- p.adjust(pval, method = "BH")
  rank <- rank(pval)
  res.ee <- as.data.frame(cbind(FDR, rank))
  rownames(res.ee) <- rownames(data)
  degs.ee <- rownames(data[FDR < param_FDR, ])
  list[["E-E(edgeR)"]][ , p] <- rownames(res.ee[order(res.ee$rank), ])

  ##  edgeR_robust
  er <- DGEList(counts = data, group = group)
  er <- edgeR::calcNormFactors(er, method = "TMM")
  er <- edgeR::estimateGLMRobustDisp(er, design = design, prior.df = 10, maxit = 6, record = FALSE)
  er <- edgeR::glmFit(er, design = design)
  er <- edgeR::glmLRT(er, coef = 2:3, test = "chisq")
  res <- er$table
  pval <- res$PValue
  pval[is.na(pval)] <- 1
  FDR <- p.adjust(pval, method = "BH")
  rank <- rank(pval)
  res.er <- as.data.frame(cbind(FDR, rank))
  rownames(res.er) <- rownames(data) 
  degs.er <- rownames(data[FDR < param_FDR, ])
  list[["edgeR_robust"]][ , p] <- rownames(res.er[order(res.er$rank), ])

  ##  D-D(DESeq)
  #dd <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
  #DE <- estimateDE(dd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
  dd <- newCountDataSet(data, group)
  dd <- DESeq::estimateSizeFactors(dd, locfunc = median)
  #sizeFactors(dd) <- sizeFactors(dd)/mean(sizeFactors(dd))
  dd <- DESeq::estimateDispersions(dd, method = "pooled", sharingMode = "maximum", fitType = "parametric")
  reduced <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ 1)
  full <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ condition)
  pval <- DESeq::nbinomGLMTest(resFull = full, resReduced = reduced)
  pval[is.na(pval)] <- 1
  FDR <- p.adjust(pval, method = "BH")
  rank <- rank(pval)
  res.dd <- as.data.frame(cbind(FDR, rank))
  rownames(res.dd) <- rownames(data) 
  degs.dd <- rownames(data[FDR < param_FDR, ])
  list[["D-D(DESeq)"]][ , p] <- rownames(res.dd[order(res.dd$rank), ])

  ##  S-S(DESeq2)
  #ss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = NULL, iteration = FALSE)
  #DE <- estimateDE(ss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
  colData <- data.frame(condition = group)
  ss <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
  ss <- DESeq2::estimateSizeFactors(ss, type = "ratio")
  #sizeFactors(ss) <- sizeFactors(ss)/mean(sizeFactors(ss))
  ss <- DESeq2::estimateDispersions(ss, fitType = "parametric")
  ss <- DESeq2::nbinomLRT(ss, full = ~ condition, reduced = ~ 1, modelMatrixType = "standard")
  res <- as.data.frame(results(ss))
  pval <- res$pvalue
  pval[is.na(pval)] <- 1
  FDR <- res$padj
  FDR[is.na(FDR)] <- 1
  rank <- rank(pval)
  res.ss <- as.data.frame(cbind(FDR, rank))
  rownames(res.ss) <- rownames(data)
  degs.ss <- rownames(data[FDR < param_FDR, ])
  list[["S-S(DESeq2)"]][ , p] <- rownames(res.ss[order(res.ss$rank), ])

  ##  voom in limma
  li <- DGEList(count = data) 
  li <- edgeR::calcNormFactors(li, method = "TMM") 
  li <- voom(li, design, plot = FALSE)
  li <- lmFit(li, design = design)    
  li <- eBayes(li)
  res <- topTable(li, coef = 2:3, n = nrow(data), sort.by = "none")
  pval <- res$P.Value
  pval[is.na(pval)] <- 1
  FDR <- res$adj.P.Val
  FDR[is.na(FDR)] <- 1
  rank <- rank(pval)
  res.li <- as.data.frame(cbind(FDR, rank))
  rownames(res.li) <- rownames(data)
  degs.li <- rownames(data[FDR < param_FDR, ])
  list[["voom"]][ , p] <- rownames(res.li[order(res.li$rank), ])

  ##  SAMseq in samr
  set.seed(2015)
  sa <- SAMseq(data, group, nperms = 100, nresamp = 20, 
               geneid = rownames(data), genenames = rownames(data),
               resp.type = "Multiclass", fdr.output = 1.0)
  sa <- rbind(sa$siggenes.table$genes.up, sa$siggenes.table$genes.lo)
  FDR <- rep(1, nrow(data))
  FDR[match(sa[,1], rownames(data))] <- as.numeric(sa[,7])/100
  rank <- rank(FDR)
  res.sa <- as.data.frame(cbind(FDR, rank))
  rownames(res.sa) <- rownames(data)
  degs.sa <- rownames(data[FDR < param_FDR, ])
  list[["SAMseq"]][ , p] <- rownames(res.sa[order(res.sa$rank), ])

  ##  PoissonSeq
  set.seed(2015)
  dat <- list(n = data, y = as.numeric(group), type = "multiclass", pair = FALSE, gname = rownames(data))
  para <- list(trans = FALSE, npermu = 500, ct.sum = -1, ct.mean = -1) 
  res <- PS.Main(dat = dat, para = para)
  res <- res[rownames(data),]
  pval <- res$pval
  pval[is.na(pval)] <- 1
  FDR <- res$fdr
  rank <- rank(pval)
  res.po <- as.data.frame(cbind(FDR, rank))
  rownames(res.po) <- rownames(data)
  degs.po <- rownames(data[FDR < param_FDR, ])
  list[["PoissonSeq"]][ , p] <- rownames(res.po[order(res.po$rank), ])

  ##  baySeq
  NDE <- factor(c(rep("NDE", Nrep * 3)))
  DE  <- factor(c(rep("G1", Nrep), rep("G2", Nrep), rep("G3", Nrep)))
  ba <- new("countData", data = data, replicates = group, groups = list(NDE = NDE, DE = DE))
  libsizes(ba) <- getLibsizes(ba, estimationType = "edgeR")
  ba <- getPriors.NB(ba, samplesize = 5000, estimation = "QL", cl = NULL)
  ba <- getLikelihoods(ba, pET = "BIC", nullData = FALSE, cl = NULL)
  res <- topCounts(ba, group = "DE", normaliseData = T, number = nrow(data))
  res <- res[rownames(data),]
  FDR <- res$FDR.DE
  rank <- rank(FDR)
  res.ba <- as.data.frame(cbind(FDR, rank))
  rownames(res.ba) <- rownames(data)
  degs.ba <- rownames(data[FDR < param_FDR, ])
  list[["baySeq"]][ , p] <- rownames(res.ba[order(res.ba$rank), ])

  ## EBSeq
  Conditions <- c(rep("G1", Nrep), rep("G2", Nrep), rep("G3", Nrep))
  PosParti <- GetPatterns(Conditions)
  MultiSize <- MedianNorm(data)
  eb <- EBMultiTest(data, NgVector = NULL, Conditions = Conditions, AllParti = PosParti,
  sizeFactors = MultiSize, maxround = 5, Qtrm = 1.0, QtrmCut = -1)
  eb <- GetMultiPP(eb)
  FDR <- eb$PP[, 1]
  rank <- rank(FDR)
  res.eb <- as.data.frame(cbind(FDR, rank))
  rownames(res.eb) <- rownames(data)
  degs.eb <- rownames(data[FDR < param_FDR, ])
  list[["EBSeq"]][ , p] <- rownames(res.eb[order(res.eb$rank), ])
}


##---------------------------------------------------------------------------##
## Sample clustering of ranked gene lists (Figure1)                          ##
##---------------------------------------------------------------------------##
data_rank <- cbind(res.eeee$rank, res.dddd$rank, res.ssss$rank, res.ee$rank, 
                   res.er$rank, res.dd$rank, res.ss$rank, res.li$rank, 
                   res.sa$rank, res.po$rank, res.ba$rank, res.eb$rank)
rownames(data_rank) <- rownames(data.ori)
colnames(data_rank) <- methods

data.dist <- as.dist(1 - cor(data_rank, method="spearman"))
hc <- hclust(data.dist, method="average")

png(file = "Figure1.png", width = 800, height = 500)
par(mar=c(0, 4, 0, 0))
plot(hc, sub="", xlab ="", ylab="Height", main="", cex=1.3, cex.lab=1.2)
dev.off()


##---------------------------------------------------------------------------##
## Numbers of DEGs and overlaps between all pairs of pipelines (Additional8) ##
## common genes                                                              ##
##---------------------------------------------------------------------------##
degs <- list(degs.eeee, degs.dddd, degs.ssss, degs.ee, 
              degs.er, degs.dd, degs.ss, degs.li,
              degs.sa, degs.po, degs.ba, degs.eb)
names(degs) <- methods
shared <- matrix(0, length(methods), length(methods))
colnames(shared) <- rownames(shared) <- methods
for(m in 1:length(methods)){
  hoge <- NULL
  for(n in 1:length(methods)){
      x <- length(intersect(degs[[m]], degs[[n]]))
      hoge <- c(hoge, x)
  }
  shared[m, ] <- hoge
}
tmp <- cbind(rownames(shared), shared)
write.table(tmp, "Additional8_Sheet1.txt", sep="\t", append=F, quote=F, row.names=F)


##---------------------------------------------------------------------------##
## Numbers of DEGs and overlaps between all pairs of pipelines (Additional8) ##
## Jaccard coefficient                                                       ##
##---------------------------------------------------------------------------##
shared <- matrix(0, length(methods), length(methods))
colnames(shared) <- rownames(shared) <- methods
for(m in 1:length(methods)){
  hoge <- NULL
  for(n in 1:length(methods)){
      x <- length(intersect(degs[[m]], degs[[n]])) / length(union(degs[[m]], degs[[n]]))
      hoge <- c(hoge, x)
  }
  shared[m, ] <- hoge
}
tmp <- cbind(rownames(shared), shared)
write.table(tmp, "Additional8_Sheet2.txt", sep="\t", append=F, quote=F, row.names=F)


##---------------------------------------------------------------------------##
## Classification of expression patterns of DEGs based on baySeq (Table4)    ##
##---------------------------------------------------------------------------##
hoge <- table(unlist(degs))
common <- names(hoge)[hoge == length(methods)]
classification <- rbind(table(res$ordering), table(res[common, ]$ordering),
                        table(res[degs.eeee, ]$ordering), table(res[degs.dddd, ]$ordering),
                        table(res[degs.ssss, ]$ordering), table(res[degs.ee, ]$ordering),
                        table(res[degs.er, ]$ordering), table(res[degs.dd, ]$ordering),
                        table(res[degs.ss, ]$ordering), table(res[degs.li, ]$ordering),
                        table(res[degs.sa, ]$ordering), table(res[degs.po, ]$ordering),
                        table(res[degs.ba, ]$ordering), table(res[degs.eb, ]$ordering))
rownames(classification) <- c("all_genes", "common", methods)
total <- rowSums(classification)
classification <- classification / total
classification <- cbind(classification, total)
tmp <- cbind(rownames(classification), classification)
write.table(tmp, "Table4.txt", sep="\t", append=F, quote=F, row.names=F)


##---------------------------------------------------------------------------##
## Classification of expression patterns of DEGs based on EBSeq (Additional9)##
##---------------------------------------------------------------------------##
hoge <- table(unlist(degs))
common <- names(hoge)[hoge == length(methods)]
classification <- rbind(table(eb$MAP), c(0, table(eb$MAP[common])),
                        table(eb$MAP[degs.eeee]), table(eb$MAP[degs.dddd]),
                        table(eb$MAP[degs.ssss]), table(eb$MAP[degs.ee]),
                        table(eb$MAP[degs.er]), table(eb$MAP[degs.dd]),
                        table(eb$MAP[degs.ss]), table(eb$MAP[degs.li]),
                        table(eb$MAP[degs.sa]), table(eb$MAP[degs.po]),
                        table(eb$MAP[degs.ba]), c(0, table(eb$MAP[degs.eb])))
rownames(classification) <- c("all_genes", "common", methods)
total <- rowSums(classification)
classification <- classification / total
classification <- cbind(classification, total)
tmp <- cbind(rownames(classification), classification)
write.table(tmp, "Additional9.txt", sep="\t", append=F, quote=F, row.names=F)

##---------------------------------------------------------------------------##
## Reproducibility for ranked gene lists (Figure2)                           ##
##---------------------------------------------------------------------------##
methods_tmp <- c("EEE-E", "DDD-D", "SSS-S", "E-E", "edgeR_robust", "D-D", 
             "S-S", "voom", "SAMseq", "PoissonSeq", "baySeq", "EBSeq")
top_num <- 100
color <- c("black", "gray", "blue", "red")
hoge <- matrix(0, nrow = 4, ncol = length(methods))
for (i in 1:length(methods)) {
  x <- list[[methods[i]]]
  hoge[1,i] <- length(intersect(head(x[,4], n=top_num), head(x[,1], n=top_num)))
  hoge[2,i] <- length(intersect(head(x[,4], n=top_num), head(x[,2], n=top_num)))
  hoge[3,i] <- length(intersect(head(x[,4], n=top_num), head(x[,3], n=top_num)))
  hoge[4,i] <- length(intersect(head(x[,1], n=top_num), intersect(head(x[,2], n=top_num), head(x[,3], n=top_num))))
}
colnames(hoge) <- methods_tmp
png("Figure2a.png", width = 1000, height = 550)
barplot(hoge, beside=TRUE, las=1, main=paste("Top", top_num, sep=""), cex.main=1.5, cex.sub=1.2,
        ylab="Number of common genes", ylim=c(0, top_num), cex.lab=1.3, col=color)
legend("topright", legend = c("rep1-6 vs. rep1-2", "rep1-6 vs. rep3-4",
 "rep1-6 vs. rep5-6", "rep1-2 vs. rep3-4 vs. rep5-6" ), col = color, pch = 15, cex = 1.2)
  abline(h=c(20, 40, 60, 80, 100), col = "gray", lty=3)
dev.off()

top_num <- 1000
color <- c("black", "gray", "blue", "red")
hoge <- matrix(0, nrow = 4, ncol = length(methods))
for (i in 1:length(methods)) {
  x <- list[[methods[i]]]
  hoge[1,i] <- length(intersect(head(x[,4], n=top_num), head(x[,1], n=top_num)))
  hoge[2,i] <- length(intersect(head(x[,4], n=top_num), head(x[,2], n=top_num)))
  hoge[3,i] <- length(intersect(head(x[,4], n=top_num), head(x[,3], n=top_num)))
  hoge[4,i] <- length(intersect(head(x[,1], n=top_num), intersect(head(x[,2], n=top_num), head(x[,3], n=top_num))))
}
colnames(hoge) <- methods_tmp
png("Figure2b.png", width = 1000, height = 550)
barplot(hoge, beside=TRUE, las=1, main=paste("Top", top_num, sep=""), cex.main=1.5, cex.sub=1.2,
        ylab="Number of common genes", ylim=c(0, top_num), cex.lab=1.3, col=color)
legend("topright", legend = c("rep1-6 vs. rep1-2", "rep1-6 vs. rep3-4",
 "rep1-6 vs. rep5-6", "rep1-2 vs. rep3-4 vs. rep5-6" ), col = color, pch = 15, cex = 1.2)
  abline(h=c(200, 400, 600, 800, 1000), col = "gray", lty=3)
dev.off()


##---------------------------------------------------------------------------##
## POG figure (Additional10)                                                 ##
##---------------------------------------------------------------------------##
size <- 500
color <- c("black", "gray", "blue", "red")
hoge <- c(1:200, seq(202, 1000, by=2), seq(1010, 10000, by=10), seq(10030, Ngene, by=30))

for (m in 1:length(methods)) {
  x <- list[[methods[m]]]
  png(file = paste("Additional10_", methods[m], ".png", sep=""), width = size, height = size)
  for(p in 1:3){
    j <- NULL
    for (i in hoge) {
      temp <- length(intersect(head(x[, 4], n = i), head(x[, p], n = i)))
      percent <- temp / i * 100
      j <- c(j, percent)
    }
    par(mar=c(4.2, 4.0, 1.2, 0.1))
    plot(hoge, j, log = "x", col = color[p], cex = 0.6, xlim = c(1, Ngene), 
      ylim = c(0, 100), ylab = '', xlab = '', axes = T)
    par(new=T) 
  }
  j <- NULL
  for (i in hoge) {
    temp <- length(intersect(head(x[, 1], n = i), intersect(head(x[, 2], n = i), head(x[, 3], n = i))))
    percent <- temp / i * 100
    j <- c(j, percent)
  }
  par(mar=c(4.2, 4.0, 1.2, 0.1))
  plot(hoge, j, log = "x", col = color[4], main = methods[m], 
    ylab = "POG(%)", xlab = "Number of top-ranked genes", cex = 0.6, xlim = c(1, Ngene), 
    ylim = c(0, 100), axes = FALSE, cex.lab=1.3, cex.main=1.3)
  axis(1, at = c(1, 10, 100, 1000, 10000))
  axis(2, at = c(0, 20, 40, 60, 80, 100))
  legend("bottomright", legend = c("rep1-6 vs. rep1-2", "rep1-6 vs. rep3-4",
 "rep1-6 vs. rep5-6", "rep1-2 vs. rep3-4 vs. rep5-6" ), col = color, pch = 1, cex = 1.2)
  abline(v = c(1, 10, 100, 1000, 10000), h = c(20, 40, 60, 80, 100), col = "gray", lty=3)
  box()
  dev.off()
}




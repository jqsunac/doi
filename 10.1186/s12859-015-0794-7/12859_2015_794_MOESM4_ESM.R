## R code for obtaining simulation results without replicates (Table  3 ).
## After execution of this R-code with default parameter settings, the
## average AUC values of 100 trials under the following conditions shown
## in Table 3 can be obtained: PDEG = 25 %, (0.5, 0.4, 0.1) for (PG1, PG2, PG3),
## and Nrep = 1. 


library(TCC) # v1.7.15

##---------------------------------------------------------------------------##
## parameter                                                                 ##
##---------------------------------------------------------------------------##
PDEG <- 0.25                    # 0.05 or 0.25
DEG.assign <- c(0.5, 0.4, 0.1)  # (PG1, PG2, PG3)
Nrep <- 1                       # Number of replicates per group
SEED <- 100                     # Number of simulation trials
methods <- c("EEE-E", "DED-E", "EDE-E", "DDD-E", "EEE-D", "DED-D", "EDE-D", "DDD-D", 
            "E-E(edgeR)", "D-E", "E-D", "D-D(DESeq)", "SSS-S", "EEE-S", "DED-S", 
            "EDE-S", "DDD-S", "S-S(DESeq2)", "E-S", "D-S")

group <- gl(3, Nrep)
design <- model.matrix(~ group)	
auc.res <- matrix(0, length(methods), SEED)
rownames(auc.res) <- methods
colnames(auc.res) <- 1:SEED

##---------------------------------------------------------------------------##
## Main process                                                              ##
##---------------------------------------------------------------------------##
for (s in 1:SEED) {
    if (s%%1 == 0) { print(paste(s, date(), sep = " ")) }
    set.seed(s)
    tcc <- simulateReadCounts(PDEG = PDEG, DEG.assign = DEG.assign, 
	DEG.foldchange = c(4, 4, 4), replicates = c(rep(Nrep, 3)))
    trueDEG <- tcc$simulation$trueDEG
    data <- tcc$count

    ## EEE-E
    eeee <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
    DE <- estimateDE(eeee, test.method = "edger",  design = design, coef = 2:3)
    auc.res[1, s] <- calcAUCValue(DE)

    ## DED-E
    dede <- calcNormFactors(tcc, norm.method = "deseq", test.method = "edger", iteration = 3)
    DE <- estimateDE(dede, test.method = "edger",  design = design, coef = 2:3)
    auc.res[2, s] <- calcAUCValue(DE)

    ## EDE-E
    edee <- calcNormFactors(tcc, norm.method = "tmm", test.method = "deseq", iteration = 3)
    DE <- estimateDE(edee, test.method = "edger",  design = design, coef = 2:3)
    auc.res[3, s] <- calcAUCValue(DE)

    ## DDD-E
    ddde <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", iteration = 3)
    DE <- estimateDE(ddde, test.method = "edger",  design = design, coef = 2:3)
    auc.res[4, s] <- calcAUCValue(DE)

    ## EEE-D
    eeed <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
    DE <- estimateDE(eeed, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[5, s] <- calcAUCValue(DE)

    ## DED-D
    dedd <- calcNormFactors(tcc, norm.method = "deseq", test.method = "edger", iteration = 3)
    DE <- estimateDE(dedd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[6, s] <- calcAUCValue(DE)

    ## EDE-D
    eded <- calcNormFactors(tcc, norm.method = "tmm", test.method = "deseq", iteration = 3)
    DE <- estimateDE(eded, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[7, s] <- calcAUCValue(DE)

    ## DDD-D
    dddd <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", iteration = 3)
    DE <- estimateDE(dddd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[8, s] <- calcAUCValue(DE)

    ## E-E (edgeR)
    #ee <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(ee, test.method = "edger",  design = design, coef = 2:3)
    #auc.res[9, s] <- calcAUCValue(DE)
    ee <- DGEList(counts = data, group = group)
    ee <- edgeR::calcNormFactors(ee, method = "TMM")
    ee <- edgeR::estimateGLMCommonDisp(ee, method = "deviance", robust = TRUE, subset = NULL)
    ee <- edgeR::glmFit(ee, design)
    ee <- edgeR::glmLRT(ee, coef = 2:3)
    r <- rank(ee$table$PValue)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[9, s] <- AUC(out)

    ## D-E
    de <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(de, test.method = "edger",  design = design, coef = 2:3)
    auc.res[10, s] <- calcAUCValue(DE)

    ## E-D
    ed <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(ed, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[11, s] <- calcAUCValue(DE)

    ## D-D (DESeq)
    #dd <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(dd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    #auc.res[12, s] <- calcAUCValue(DE)
    dd <- newCountDataSet(data, group)
    dd <- DESeq::estimateSizeFactors(dd, locfunc = median)
    #sizeFactors(dd) <- sizeFactors(dd)/mean(sizeFactors(dd))
    dd <- DESeq::estimateDispersions(dd, method = "blind", sharingMode = "fit-only", fitType = "parametric")
    reduced <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ 1)
    full <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ condition)
    pval <- DESeq::nbinomGLMTest(resFull = full, resReduced = reduced)
    pval[is.na(pval)] <- 1
    r <- rank(pval)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[12, s] <- AUC(out)

    ## SSS-S
    ssss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2", iteration = 3)
    DE <- estimateDE(ssss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[13, s] <- calcAUCValue(DE)

    ## EEE-S
    eees <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
    DE <- estimateDE(eees, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[14, s] <- calcAUCValue(DE)

    ## DED-S
    deds <- calcNormFactors(tcc, norm.method = "deseq", test.method = "edger", iteration = 3)
    DE <- estimateDE(deds, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[15, s] <- calcAUCValue(DE)

    ## EDE-S
    edes <- calcNormFactors(tcc, norm.method = "tmm", test.method = "deseq", iteration = 3)
    DE <- estimateDE(edes, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[16, s] <- calcAUCValue(DE)

    ## DDD-S
    ddds <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", iteration = 3)
    DE <- estimateDE(ddds, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[17, s] <- calcAUCValue(DE)

    ## S-S (DESeq2)
    #ss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(ss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    #auc.res[18, s] <- calcAUCValue(DE)
    colData <- data.frame(condition = group)
    ss <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
    ss <- DESeq2::estimateSizeFactors(ss, type = "ratio")
    #sizeFactors(ss) <- sizeFactors(ss)/mean(sizeFactors(ss))
    ss <- DESeq2::estimateDispersions(ss, fitType = "parametric")
    ss <- DESeq2::nbinomLRT(ss, full = ~ condition, reduced = ~ 1, modelMatrixType = "standard")
    res <- results(ss)
    res$pvalue[is.na(res$pvalue)] <- 1
    r <- rank(res$pvalue)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[18, s] <- AUC(out)

    ## E-S
    es <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(es, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[19, s] <- calcAUCValue(DE)

    ## D-S
    ds <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(ds, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    auc.res[20, s] <- calcAUCValue(DE)
}
##---------------------------------------------------------------------------##
## Table 3                                                                   ##
##---------------------------------------------------------------------------##
auc.ave <- rowMeans(auc.res)
auc.ave

tmp <- cbind(rownames(auc.res), auc.res)
filename <- paste("Table3", PDEG, DEG.assign[1], DEG.assign[2], "raw.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)
tmp <- cbind(names(auc.ave), auc.ave)
filename <- paste("Table3", PDEG, DEG.assign[1], DEG.assign[2], "ave.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)

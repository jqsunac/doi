## R code for obtaining simulation results with replicates (Table  1 ).
## After execution of this R-code with default parameter settings, the
## average AUC values of 100 trials under the following conditions shown
## in Table 1 can be obtained: PDEG = 5 %, (0.5, 0.4, 0.1) for (PG1, PG2, PG3),
## and Nrep = 3. The results in Additional file 2 can also be obtained by
## changing the parameter Nrep to be 6 or 9.


library(TCC)		# v1.7.15
library(samr)		# v2.0
library(PoissonSeq)	# v1.1.2
library(EBSeq)		# v1.6.0
library(limma)		# v3.22.2

##---------------------------------------------------------------------------##
## parameter                                                                 ##
##---------------------------------------------------------------------------##
PDEG <- 0.05                    # 0.05 or 0.25
#DEG.assign <- c(1/3, 1/3, 1/3)  # (PG1, PG2, PG3)
#DEG.assign <- c(0.5, 0.3, 0.2)  # (PG1, PG2, PG3)
#DEG.assign <- c(0.5, 0.4, 0.1)  # (PG1, PG2, PG3)
#DEG.assign <- c(0.6, 0.2, 0.2)  # (PG1, PG2, PG3)
#DEG.assign <- c(0.6, 0.3, 0.1)  # (PG1, PG2, PG3)
#DEG.assign <- c(0.7, 0.2, 0.1)  # (PG1, PG2, PG3)
DEG.assign <- c(0.8, 0.1, 0.1)  # (PG1, PG2, PG3)
Nrep <- 3                       # Number of replicates per group
SEED <- 2                     # Number of simulation trials
methods <- c("EEE-E(TCC)", "DDD-D(TCC)", "SSS-S(TCC)", "E-E(edgeR)", "edgeR_robust", 
         "D-D(DESeq)", "S-S(DESeq2)", "voom", "SAMseq", "PoissonSeq", "baySeq", "EBSeq") 
filename_prefix <- "result"     #

group <- gl(3, Nrep)
design <- model.matrix(~ group)	
auc.res <- matrix(0, length(methods), SEED)
rownames(auc.res) <- methods
colnames(auc.res) <- 1:SEED

time.res <- matrix(0, length(methods), SEED)
rownames(time.res) <- methods
colnames(time.res) <- 1:SEED

pauc.res <- matrix(0, length(methods), SEED)
rownames(pauc.res) <- methods
colnames(pauc.res) <- 1:SEED

##---------------------------------------------------------------------------##
## Main process                                                              ##
##---------------------------------------------------------------------------##
for (s in 1:SEED){
    if(s%%1 == 0) { print(paste(s, date(), sep=" ")) }
    set.seed(s)
    #####################################
    ### Generation of simulation data ###
    #####################################
    ### a fixed level of DE (four-fold for all DEGs)
    ### For obtaining Additinal file 1 Sheet 1-4
    tcc <- simulateReadCounts(PDEG = PDEG, DEG.assign = DEG.assign, 
            DEG.foldchange = c(4, 4, 4), replicates = c(rep(Nrep, 3)))
    trueDEG <- as.numeric(tcc$simulation$trueDEG != 0)
    data <- tcc$count

    ### different levels and distributions of DE, where
    ### the fold-changes for DEGs are randomly sampled from 
    ### "1.2 + a gamma distribution with shape = 2.0 and scale = 0.5"
    ### mean fold-change is 2.2 (= 1.2 + 2.0 * 0.5)
    ### For obtaining Additinal file 1 Sheet 5
    #source("http://www.iu.a.u-tokyo.ac.jp/~kadota/TCC/TCC.simulation.R")
    #fc.params <- data.frame(
    #                   floor = c(1.2, 1.2, 1.2),
    #                   shape = c(2.0, 2.0, 2.0),
    #                   scale = c(0.5, 0.5, 0.5))
    #fcm <- makeFCMatrix(Ngene = 10000, PDEG = PDEG, 
    #         replicates = c(rep(Nrep, 3)), DEG.assign = DEG.assign,
    #         fc.params = fc.params)
    #mean(fcm[fcm != 1])
    #tcc <- simulateReadCounts(Ngene = 10000,
    #         PDEG = PDEG, DEG.assign = DEG.assign, 
    #         DEG.foldchange = NULL, replicates = c(rep(Nrep, 3)),
    #         fc.matrix = fcm)
    #trueDEG <- as.numeric(tcc$simulation$trueDEG != 0)
    #data <- tcc$count

    ########################################
    ### Differential Expression Analysis ###
    ########################################
    ### EEE-E
    ptm <- proc.time()
    eeee <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)
    DE <- estimateDE(eeee, test.method = "edger", design = design, coef = 2:3)
    time.res[1, s] <- (proc.time() - ptm)[3]
    r <- rank(getResult(DE)$p.value)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[1, s] <- AUC(out)
    pauc.res[1, s] <- pAUC(out, t0 = 0.1)

    ### DDD-D
    ptm <- proc.time()
    dddd <- calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq", iteration = 3)
    DE <- estimateDE(dddd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    time.res[2, s] <- (proc.time() - ptm)[3]
    r <- rank(getResult(DE)$p.value)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[2, s] <- AUC(out)
    pauc.res[2, s] <- pAUC(out, t0 = 0.1)

    ### SSS-S
    ptm <- proc.time()
    ssss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2", iteration = 3)
    DE <- estimateDE(ssss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    time.res[3, s] <- (proc.time() - ptm)[3]
    r <- rank(getResult(DE)$p.value)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[3, s] <- AUC(out)
    pauc.res[3, s] <- pAUC(out, t0 = 0.1)

    ### E-E (edgeR)
    ptm <- proc.time()
    #ee <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(ee, test.method = "edger", design = design, coef = 2:3)
    #auc.res[4, s] <- calcAUCValue(DE)
    ee <- DGEList(counts = data, group = group)
    ee <- edgeR::calcNormFactors(ee, method = "TMM")
    ee <- edgeR::estimateGLMCommonDisp(ee, design, method = "CoxReid")
    ee <- edgeR::estimateGLMTrendedDisp(ee, design, method = "auto")
    ee <- edgeR::estimateGLMTagwiseDisp(ee, design)
    ee <- edgeR::glmFit(ee, design)
    ee <- edgeR::glmLRT(ee, coef = 2:3, test = "chisq")
    time.res[4, s] <- (proc.time() - ptm)[3]
    r <- rank(ee$table$PValue)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[4, s]<- AUC(out)
    pauc.res[4, s] <- pAUC(out, t0 = 0.1)

    ### edgeR_robust
    ptm <- proc.time()
    er <- DGEList(counts = data, group = group)
    er <- edgeR::calcNormFactors(er, method = "TMM")
    er <- estimateGLMRobustDisp(er, design = design, prior.df = 10, maxit = 6, record = FALSE)
    er <- edgeR::glmFit(er, design = design)
    er <- edgeR::glmLRT(er, coef = 2:3, test = "chisq")
    time.res[5, s] <- (proc.time() - ptm)[3]
    pval <- er$table$PValue
    r <- rank(pval)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[5, s] <- AUC(out)
    pauc.res[5, s] <- pAUC(out, t0 = 0.1)

    ### D-D (DESeq)
    ptm <- proc.time()
    #dd <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(dd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    #auc.res[6, s] <- calcAUCValue(DE)
    dd <- newCountDataSet(data, group)
    dd <- DESeq::estimateSizeFactors(dd, locfunc = median)
    #sizeFactors(dd) <- sizeFactors(dd)/mean(sizeFactors(dd))
    dd <- DESeq::estimateDispersions(dd, method = "pooled", sharingMode = "maximum", fitType = "parametric")
    reduced <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ 1)
    full <- DESeq::fitNbinomGLMs(dd, modelFormula = count ~ condition)
    pval <- DESeq::nbinomGLMTest(resFull = full, resReduced = reduced)
    time.res[6, s] <- (proc.time() - ptm)[3]
    pval[is.na(pval)] <- 1
    r <- rank(pval)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[6, s] <- AUC(out)
    pauc.res[6, s] <- pAUC(out, t0 = 0.1)

    ### S-S (DESeq2)
    ptm <- proc.time()
    #ss <- calcNormFactors(tcc, norm.method = "deseq2", test.method = NULL, iteration = FALSE)
    #DE <- estimateDE(ss, test.method = "deseq2", full = ~ group, reduced = ~ 1)
    #auc.res[7, s] <- calcAUCValue(DE)
    colData <- data.frame(condition = group)
    ss <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
    ss <- DESeq2::estimateSizeFactors(ss, type = "ratio")
    #sizeFactors(ss) <- sizeFactors(ss)/mean(sizeFactors(ss))
    ss <- DESeq2::estimateDispersions(ss, fitType = "parametric")
    ss <- DESeq2::nbinomLRT(ss, full = ~ condition, reduced = ~ 1, modelMatrixType = "standard")
    res <- results(ss)
    time.res[7, s] <- (proc.time() - ptm)[3]
    res$pvalue[is.na(res$pvalue)] <- 1
    r <- rank(res$pvalue)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[7, s] <- AUC(out)
    pauc.res[7, s] <- pAUC(out, t0 = 0.1)

    ### voom in limma
    ptm <- proc.time()
    li <- DGEList(count = data)
    li <- edgeR::calcNormFactors(li, method = "TMM") 
    li <- voom(li, design, plot = FALSE)
    li <- lmFit(li, design = design) 
    li <- eBayes(li)
    res <- topTable(li, coef = 2:3, n = nrow(data), sort.by = "none")
    time.res[8, s] <- (proc.time() - ptm)[3]
    r <- rank(res$P.Value)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[8, s] <- AUC(out)
    pauc.res[8, s] <- pAUC(out, t0 = 0.1)

    ### SAMseq in samr
    set.seed(2015)
    ptm <- proc.time()
    sa <- SAMseq(data, tcc$group$group, nperms = 100, nresamp = 20, 
                 geneid = rownames(data), genenames = rownames(data),
                 resp.type = "Multiclass", fdr.output = 1.0)
    sa <- rbind(sa$siggenes.table$genes.up, sa$siggenes.table$genes.lo)
    time.res[9, s] <- (proc.time() - ptm)[3]
    FDR <- rep(1, nrow(data))
    FDR[match(sa[,1], rownames(data))] <- as.numeric(sa[,7])/100
    r <- rank(FDR)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[9, s] <- AUC(out)
    pauc.res[9, s] <- pAUC(out, t0 = 0.1)

    ### PoissonSeq
    set.seed(2015)
    ptm <- proc.time()
    dat <- list(n = data, y = as.numeric(group), type = "multiclass", pair = FALSE, gname = rownames(data))
    para <- list(trans = FALSE, npermu = 500, ct.sum = -1, ct.mean = -1)
    res <- PS.Main(dat = dat, para = para)
    time.res[10, s] <- (proc.time() - ptm)[3]
    res <- res[rownames(data), ]
    res$pval[is.na(res$pval)] <- 1
    r <- rank(res$pval)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[10, s] <- AUC(out)
    pauc.res[10, s] <- pAUC(out, t0 = 0.1)

    ### baySeq
    ptm <- proc.time()
    DE <- c(rep("G1", Nrep), rep("G2", Nrep), rep("G3", Nrep))
    NDE <- factor(c(rep("NDE", Nrep*3)))
    ba <- new("countData", data = data, replicates = group, groups = list(NDE = NDE, DE = DE))
    libsizes(ba) <- getLibsizes(ba, estimationType = "edgeR") 
    ba <- getPriors.NB(ba, samplesize = 5000, estimation = "QL", cl = NULL)
    ba <- getLikelihoods(ba, pET = "BIC", nullData = FALSE, cl = NULL)
    res <- topCounts(ba, group = "DE", normaliseData = T, number = nrow(data))
    time.res[11, s] <- (proc.time() - ptm)[3]
    res <- res[rownames(data),]
    FDR <- res$FDR.DE
    r <- rank(FDR)
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[11, s] <- AUC(out)
    pauc.res[11, s] <- pAUC(out, t0 = 0.1)

    ### EBSeq
    ptm <- proc.time()
    Conditions <- c(rep("G1", Nrep), rep("G2", Nrep), rep("G3", Nrep))
    PosParti <- GetPatterns(Conditions)
    MultiSize <- MedianNorm(data)
    eb <- EBMultiTest(data, NgVector = NULL, Conditions = Conditions, AllParti = PosParti,
    sizeFactors = MultiSize, maxround = 5, Qtrm = 1.0, QtrmCut = -1)
    eb <- GetMultiPP(eb)
    time.res[12, s] <- (proc.time() - ptm)[3]
    r <- rank(eb$PP[, 1])
    out <- rocdemo.sca(truth = trueDEG, data = -r)
    auc.res[12, s]<- AUC(out)
    pauc.res[12, s] <- pAUC(out, t0 = 0.1)
}


##---------------------------------------------------------------------------##
## Output                                                                    ##
##---------------------------------------------------------------------------##
auc.ave <- rowMeans(auc.res)
auc.ave
tmp <- cbind(rownames(auc.res), auc.res)
filename <- paste("result", "AUC", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "raw.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)
tmp <- cbind(names(auc.ave), auc.ave)
filename <- paste("result", "AUC", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "ave.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)

pauc.ave <- rowMeans(pauc.res)
pauc.ave
tmp <- cbind(rownames(pauc.res), pauc.res)
filename <- paste("result", "pAUC", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "raw.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)
tmp <- cbind(names(pauc.ave), pauc.ave)
filename <- paste("result", "pAUC", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "ave.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)

time.ave <- rowMeans(time.res)
time.ave
tmp <- cbind(rownames(time.res), time.res)
filename <- paste("result", "TIME", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "raw.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)
tmp <- cbind(names(time.ave), time.ave)
filename <- paste("result", "TIME", Nrep, PDEG, DEG.assign[1], DEG.assign[2], "ave.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)

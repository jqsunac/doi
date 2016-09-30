## R code for obtaining simulation results with replicates (Table  2 ).
## After execution of this R-code with default parameter settings,
## the average AUC values of 100 trials under the following conditions
## shown in Table 2 can be obtained: PDEG = 5 %, (0.5, 0.4, 0.1) for
## (PG1, PG2, PG3), and Nrep = 3.


library(TCC) # v1.7.15

##---------------------------------------------------------------------------##
## paramter                                                                  ##
##---------------------------------------------------------------------------##
PDEG <- 0.05                    # 0.05 or 0.25
DEG.assign <- c(0.5, 0.4, 0.1)  # (PG1, PG2, PG3)
Nrep <- 3                       # Number of replicates per group
SEED <- 100                     # Number of simulation trials
method <- c("EEE-E", "DED-E", "EDE-E", "DDD-E", "EEE-D", "DED-D", "EDE-D",
            "DDD-D", "E-E (edgeR)", "D-E", "E-D", "D-D (DESeq)")

group <- gl(3, Nrep)
design <- model.matrix(~ group)	
auc.res <- matrix(0, length(method), SEED)
rownames(auc.res) <- method
colnames(auc.res) <- 1:SEED

##---------------------------------------------------------------------------##
## Main process                                                              ##
##---------------------------------------------------------------------------##
for(s in 1:SEED){
    if(s%%1 == 0) { print(paste(s, date(), sep=" ")) }
    set.seed(s)
    tcc <- simulateReadCounts(PDEG = PDEG, DEG.assign = DEG.assign, 
        DEG.foldchange = c(4, 4, 4), replicates = c(rep(Nrep, 3)))
    trueDEG <- tcc$simulation$trueDEG

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

    ## E-E
    ee <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(ee, test.method = "edger",  design = design, coef = 2:3)
    auc.res[9, s] <- calcAUCValue(DE)

    ## D-E
    de <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(de, test.method = "edger",  design = design, coef = 2:3)
    auc.res[10, s] <- calcAUCValue(DE)

    ## E-D
    ed <- calcNormFactors(tcc, norm.method = "tmm", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(ed, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[11, s] <- calcAUCValue(DE)

    ## D-D
    dd <- calcNormFactors(tcc, norm.method = "deseq", test.method = NULL, iteration = FALSE)
    DE <- estimateDE(dd, test.method = "deseq", full = count ~ condition, reduced = count ~ 1)
    auc.res[12, s] <- calcAUCValue(DE)
}
##---------------------------------------------------------------------------##
## Table 2                                                                   ##
##---------------------------------------------------------------------------##
auc.ave <- rowMeans(auc.res)
auc.ave

tmp <- cbind(rownames(auc.res), auc.res)
filename <- paste("Table2", PDEG, DEG.assign[1], DEG.assign[2], "raw.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)
tmp <- cbind(names(auc.ave), auc.ave)
filename <- paste("Table2", PDEG, DEG.assign[1], DEG.assign[2], "ave.txt", sep="_")
write.table(tmp, filename, sep="\t", append=F, quote=F, row.names=F)

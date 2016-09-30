# simulateScript.R
#
# This script is for simulation study.
# The simulation condition need describing in arguments.
# If one executes simulation 100 times with PDEG = 20% and PG1 = 80%,
# run this script as:
#   R --slave --vanilla --args 0.2 0.8 100 < simulateScript.R
#
# It needs to set the directory path for saving simulation result.
result.path <- "./result_revise"
# The results are saved with the prefix "simData_".
#
#
# Reuqirement:
#    R 2.15.2 (2012-10-26)
#    TCC 1.0.0
#    ROC 1.34.0
#    DESeq 1.10.1
#    baySeq 1.12.0
#    edgeR 3.0.3
#
#




# load packages.
library(TCC)

# arguments analysis.
#args <- commandArgs()
#PDEG <- c(as.numeric(args[5]))
#PA <- c(as.numeric(args[6]))
#TIMES <- c(as.numeric(args[7]))
PDEG <- 0.2
PA <- 0.9
TIMES <- 3

simulateMultiReps <- function(path, times, PDEG, PA) {
  # simulate study for estimation accuracy.
  # 1) Generate simulation data.
  # 2) Normalize data with TMM, TbT, DES and iDES2.
  # 3) Estimate DEGs with four types normlized count.
  # 4) Caculate AUC values of four types.
  tmpl <- matrix(0, nrow = times, ncol = 4)  # template
  colnames(tmpl) <- c("TMM", "TbT", "DEGES", "iDEGES3")
  acc <- list(TP = tmpl, TN = tmpl, FP = tmpl, FN = tmpl)
  time <- tmpl
  auc <- tmpl
  threshold.type <- tmpl

  arg.des  <- list(FALSE,   TRUE,     TRUE,    3)
  arg.norm <- c("tmm",   "tmm",    "tmm",   "tmm")
  arg.test <- c("edger", "bayseq", "edger", "edger")

  # Sample simulation data.
  for (i in 1:times) {
    tcc <- generateSimulationData(PDEG = PDEG, DEG.assign = c(PA, 1 - PA))
    for (j in 1:length(arg.norm)) {
      tcc <- calcNormFactors(tcc, norm.method = arg.norm[j],
        test.method = arg.test[j], iteration = arg.des[[j]], samplesize = 10)
      threshold.type[i, j] <- tcc$private$DEGES.threshold.type
      tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
      auc[i, j] <- calcAUCValue(tcc)
      time[i, j] <- tcc$stat$execution.time[3]
      acc["TP"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG != 0) & 
        (tcc$private$DEGES.potentialDEG != 0))) / length(tcc$simulation$trueDEG)
      acc["TN"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG == 0) & 
        (tcc$private$DEGES.potentialDEG == 0))) / length(tcc$simulation$trueDEG)
      acc["FP"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG == 0) &
        (tcc$private$DEGES.potentialDEG != 0))) / length(tcc$simulation$trueDEG)
      acc["FN"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG != 0) & 
        (tcc$private$DEGES.potentialDEG == 0))) / length(tcc$simulation$trueDEG)
      write.table(time,
        file = paste(path, "/simData_multirep_time_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep = ""),
        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      write.table(auc,
        file = paste(path, "/simData_multirep_auc_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep = ""),
        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      write.table(threshold.type,
        file = paste(path, "/simData_multirep_thresholdType_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep = ""),
        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      for (k in 1:length(acc)) {
        write.table(acc[[k]],
          file = paste(path, "/simData_multirep_", names(acc[k]), "_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep = ""),
          col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      }
    }
    remove(tcc)
  }
}


simulateMonoRep <- function(path, times, PDEG, PA) {
  # simulate study for estimation accuracy.
  # for no replicate samples.
  # 1) Generate simulation data.
  # 2) Normalize data with DESeq, DES and iDES2.
  # 3) Estimate DEGs with three types normlized count.
  # 4) Caculate AUC values of three types.
  tmpl <- matrix(0, nrow=times, ncol=3)  # template
  colnames(tmpl) <- c("DESeq", "DEGES", "iDEGES3")
  acc <- list(TP = tmpl, TN = tmpl, FP = tmpl, FN = tmpl)
  time <- tmpl
  auc <- tmpl

  arg.des  <- list(FALSE,  TRUE,  3)
  arg.norm <- c("deseq", "deseq", "deseq")
  arg.test <- c("deseq", "deseq", "deseq")

  # Sample simulation data.
  for (i in 1:times) {
    tcc <- generateSimulationData(PDEG=PDEG, DEG.assign=c(PA, 1 - PA), group=c(1,1))
    for (j in 1:length(arg.norm)) {
      tcc <- calcNormFactors(tcc,
        norm.method=arg.norm[j],
        test.method=arg.test[j],
        iteration=arg.des[[j]],
        floorPDEG=0.05)
      tcc <- estimateDE(tcc, test.method="deseq", FDR=0.1)
      auc[i, j] <- calcAUCValue(tcc)
      time[i, j] <- tcc$stat$execution.time[3]
      acc["TP"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG != 0) & 
        (tcc$private$DEGES.potentialDEG != 0))) / length(tcc$simulation$trueDEG)
      acc["TN"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG == 0) & 
        (tcc$private$DEGES.potentialDEG == 0))) / length(tcc$simulation$trueDEG)
      acc["FP"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG == 0) &
        (tcc$private$DEGES.potentialDEG != 0))) / length(tcc$simulation$trueDEG)
      acc["FN"][[1]][i, j] <- sum(as.numeric((tcc$simulation$trueDEG != 0) & 
        (tcc$private$DEGES.potentialDEG == 0))) / length(tcc$simulation$trueDEG)
    }
    remove(tcc)
  }

  write.table(time,
    file=paste(path, "/simData_monorep_time_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep=""),
    col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  write.table(auc,
    file=paste(path, "/simData_monorep_auc_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep=""),
    col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  for (i in 1:length(acc)) {
    write.table(acc[[i]],
      file=paste(path, "/simData_monorep_", names(acc[i]), "_PDEG_", PDEG, "_PG1_", PA, "_raw.txt", sep=""),
      col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  }
}

if (0) {
checkCovergence <- function(result.path, times.sim, PDEG, PA, times.try) {
  times.try <- times.try + 1
  for (t in 1:times.sim) {
    auc <- matrix(0, nrow = times.try, ncol = 2)
    colnames(auc) <- c("iDEGES/edgeR", "iDEGES/TbT")
    rownames(auc) <- 1:times.try
    
    for (i in 1:times.try) {
      tcc <- generateSimulationData(PDEG = PDEG, DEG.assign = c(PA, 1 - PA))
      tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = i - 1,
                                 FDR = 0.1, floorPDEG = 0.05)
      tcc <- estimateDE(tcc, test.method = "edger")
      auc[i, 1] <- calcAUCValue(tcc)
      tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq", iteration = i - 1,
                                 samplesize = 10000)
      tcc <- estimateDE(tcc, test.method = "edger")
      auc[i, 2] <- calcAUCValue(tcc)
      write.table(auc,
        file = paste(result.path, "/simData_multirep_limit_PDEG_", PDEG, "_PG1_", 
          PA, "_TIMES_", t, ".txt", sep = ""), 
        quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE)
    }
  }
}
}

#simulateMultiReps(result.path, TIMES, PDEG, PA)
#simulateMonoRep(result.path, TIMES, PDEG, PA)




crossValid <- function(PDEG = 0.2, PA = 0.9, cycles = 2) {
  norm.norm <- c("edger", "deseq", "tmm", "edger", "edger", "edger", "tmm", "deseq")
  norm.test <- c("edger", "deseq", "bayseq", "edger", "edger", "edger", "bayseq", "deseq")
  norm.iter <- c(0, 0, 0, 1, 1, 1, 1, 1)
  test.test <- c("edger", "deseq", "bayseq", "edger", "bayseq", "deseq", "bayseq", "deseq")

  auc <- matrix(0, ncol = length(norm.norm), nrow = cycles)
  colnames(auc) <- c("edgeR", "DESeq", "baySeq", "DEGES/edgeR-edgeR",
                     "DEGSE/edgeR-baySeq", "DEGES/edgeR-DESeq",
                     "DEGES/TbT-baySeq", "DEGES/DESeq-DESeq")

  for (n in 1:cycles) {
    tcc <- generateSimulationData(PDEG = PDEG, DEG.assign = c(PA, 1 - PA))
    for (i in 1:length(norm.norm)) {
      t <- tcc$copy()
      t <- calcNormFactors(t, norm.method = norm.norm[i],
                           test.method = norm.test[i], iteration = norm.iter[i],
                           samplesize = 10000)
      t <- estimateDE(t, test.method = test.test[i], samplesize = 10000)
      auc[n, i] <- calcAUCValue(t)
      write.table(auc, file = paste(result.path, "/simData_multirep_CrossValid",
                  "_PDEG_", PDEG, "_PG1_", PA, ".txt", sep = ""),
                  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    }
  }
}
#crossValid(PDEG = 0.2, PA = 1, cycles = 50)
#crossValid(PDEG = 0.2, PA = 0.9, cycles = 50)
#crossValid(PDEG = 0.2, PA = 0.8, cycles = 50)
#crossValid(PDEG = 0.2, PA = 0.7, cycles = 50)
#crossValid(PDEG = 0.2, PA = 0.6, cycles = 50)
#crossValid(PDEG = 0.2, PA = 0.5, cycles = 50)



crossValidKadota <- function(PDEG = 0.2, PA = 0.9, cycles = 2, iteration = 2, n = 0) {
  #cycles <- 1
  #iteration <- 8
  #load("../realdata/refData/rdr6_wt.RData")
  #d <- as.matrix(cD@data)
  #rownames(d) <- c(1:nrow(cD@data))
  #tmpl <- matrix(0, nrow = nrow(d), ncol = iteration + 2)

  tmpl <- matrix(0, nrow = 10000, ncol = iteration + 2)
  colnames(tmpl) <- c(1:(iteration + 1), "trueDEG")

  #for (n in 1:cycles) {
    potentialDEG.edger <- tmpl
    potentialDEG.tbt <- tmpl
    potentialDEG.deseq <- tmpl
    potentialDEG.tdt <- tmpl
    potentialRank.edger <- tmpl
    potentialRank.tbt <- tmpl
    potentialRank.deseq <- tmpl
    potentialRank.tdt <- tmpl
    tcc <- generateSimulationData(Ngene = 10000, PDEG = PDEG, DEG.assign = c(PA, 1 - PA))

    tcc.edger <- tcc$copy()
    tcc.tbt   <- tcc$copy()
    tcc.deseq <- tcc$copy()
    tcc.tdt   <- tcc$copy()
    # DEGSE/edgeR
    tcc.edger <- calcNormFactors(tcc.edger, norm.method = "tmm", test.method = "edger", 
                                 iteration = iteration)
    tcc.edger <- estimateDE(tcc.edger, test.method = "edger")
    for (i in 1:iteration) {
      potentialDEG.edger[, i] <- tcc.edger$private$debug$potentialDEG[[i]]
      potentialRank.edger[, i] <- tcc.edger$private$debug$potentialRank[[i]]
    }
    potentialDEG.edger[, iteration + 1] <- tcc.edger$estimatedDEG
    potentialRank.edger[, iteration + 1] <- tcc.edger$stat$rank
    potentialDEG.edger[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.edger[, iteration + 2] <- tcc$simulation$trueDEG

    write.table(potentialDEG.edger, file = paste(result.path, "/simResult_multi_Kdt_DEGESedgeR_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.edger, file = paste(result.path, "/simResult_multi_Kdt_DEGESedgeR_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    # DEGES/TbT
    tcc.tbt <- calcNormFactors(tcc.tbt, norm.method = "tmm", test.method = "bayseq", 
                                 iteration = iteration, samplesize = 10000)
    tcc.tbt <- estimateDE(tcc.tbt, test.method = "edger")
    for (i in 1:iteration) {
      potentialDEG.tbt[, i] <- tcc.tbt$private$debug$potentialDEG[[i]]
      potentialRank.tbt[, i] <- tcc.tbt$private$debug$potentialRank[[i]]
    }
    potentialDEG.tbt[, iteration + 1] <- tcc.tbt$estimatedDEG
    potentialRank.tbt[, iteration + 1] <- tcc.tbt$stat$rank
    potentialDEG.tbt[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.tbt[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.tbt, file = paste(result.path, "/simResult_multi_Kdt_DEGESTbT_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.tbt, file = paste(result.path, "/simResult_multi_Kdt_DEGESTbT_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    # DEGES/DESeq
    tcc.deseq <- calcNormFactors(tcc.deseq, norm.method = "deseq", test.method = "deseq", 
                                 iteration = iteration)
    tcc.deseq <- estimateDE(tcc.deseq, test.method = "edger")
    for (i in 1:iteration) {
      potentialDEG.deseq[, i] <- tcc.deseq$private$debug$potentialDEG[[i]]
      potentialRank.deseq[, i] <- tcc.deseq$private$debug$potentialRank[[i]]
    }
    potentialDEG.deseq[, iteration + 1] <- tcc.deseq$estimatedDEG
    potentialRank.deseq[, iteration + 1] <- tcc.deseq$stat$rank
    potentialDEG.deseq[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.deseq[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.deseq, file = paste(result.path, "/simResult_multi_Kdt_DEGESDESeq_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.deseq, file = paste(result.path, "/simResult_multi_Kdt_DEGESDESeq_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    # DEGES/TDT
    tcc.tdt <- calcNormFactors(tcc.tdt, norm.method = "tmm", test.method = "deseq", 
                                 iteration = iteration)
    tcc.tdt <- estimateDE(tcc.tdt, test.method = "edger")
    for (i in 1:iteration) {
      potentialDEG.tdt[, i] <- tcc.tdt$private$debug$potentialDEG[[i]]
      potentialRank.tdt[, i] <- tcc.tdt$private$debug$potentialRank[[i]]
    }
    potentialDEG.tdt[, iteration + 1] <- tcc.tdt$estimatedDEG
    potentialRank.tdt[, iteration + 1] <- tcc.tdt$stat$rank
    potentialDEG.tdt[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.tdt[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.tdt, file = paste(result.path, "/simResult_multi_Kdt_DEGESTDT_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.tdt, file = paste(result.path, "/simResult_multi_Kdt_DEGESTDT_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)


    potentialDEG.tbt <- tmpl
    potentialDEG.deseq <- tmpl
    potentialDEG.tdt <- tmpl
    potentialRank.tbt <- tmpl
    potentialRank.deseq <- tmpl
    potentialRank.tdt <- tmpl
    tcc <- generateSimulationData(Ngene = 10000, PDEG = PDEG, DEG.assign = c(PA, 1 - PA), group = c(1, 2))

    tcc.tbt   <- tcc$copy()
    tcc.deseq <- tcc$copy()
    tcc.tdt   <- tcc$copy()
    # DEGES/TbT
    tcc.tbt <- calcNormFactors(tcc.tbt, norm.method = "tmm", test.method = "bayseq", 
                                 iteration = iteration, samplesize = 10000)
    tcc.tbt <- estimateDE(tcc.tbt, test.method = "deseq")
    for (i in 1:iteration) {
      potentialDEG.tbt[, i] <- tcc.tbt$private$debug$potentialDEG[[i]]
      potentialRank.tbt[, i] <- tcc.tbt$private$debug$potentialRank[[i]]
    }
    potentialDEG.tbt[, iteration + 1] <- tcc.tbt$estimatedDEG
    potentialRank.tbt[, iteration + 1] <- tcc.tbt$stat$rank
    potentialDEG.tbt[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.tbt[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.tbt, file = paste(result.path, "/simResult_mono_Kdt_DEGESTbT_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.tbt, file = paste(result.path, "/simResult_mono_Kdt_DEGESTbT_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    # DEGES/DESeq
    tcc.deseq <- calcNormFactors(tcc.deseq, norm.method = "deseq", test.method = "deseq", 
                                 iteration = iteration)
    tcc.deseq <- estimateDE(tcc.deseq, test.method = "deseq")
    for (i in 1:iteration) {
      potentialDEG.deseq[, i] <- tcc.deseq$private$debug$potentialDEG[[i]]
      potentialRank.deseq[, i] <- tcc.deseq$private$debug$potentialRank[[i]]
    }
    potentialDEG.deseq[, iteration + 1] <- tcc.deseq$estimatedDEG
    potentialRank.deseq[, iteration + 1] <- tcc.deseq$stat$rank
    potentialDEG.deseq[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.deseq[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.deseq, file = paste(result.path, "/simResult_mono_Kdt_DEGESDESeq_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.deseq, file = paste(result.path, "/simResult_mono_Kdt_DEGESDESeq_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    # DEGES/TDT
    tcc.tdt <- calcNormFactors(tcc.tdt, norm.method = "tmm", test.method = "deseq", 
                                 iteration = iteration)
    tcc.tdt <- estimateDE(tcc.tdt, test.method = "deseq")
    for (i in 1:iteration) {
      potentialDEG.tdt[, i] <- tcc.tdt$private$debug$potentialDEG[[i]]
      potentialRank.tdt[, i] <- tcc.tdt$private$debug$potentialRank[[i]]
    }
    potentialDEG.tdt[, iteration + 1] <- tcc.tdt$estimatedDEG
    potentialRank.tdt[, iteration + 1] <- tcc.tdt$stat$rank
    potentialDEG.tdt[, iteration + 2] <- tcc$simulation$trueDEG
    potentialRank.tdt[, iteration + 2] <- tcc$simulation$trueDEG
    write.table(potentialDEG.tdt, file = paste(result.path, "/simResult_mono_Kdt_DEGESTDT_pDEG_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)
    write.table(potentialRank.tdt, file = paste(result.path, "/simResult_mono_Kdt_DEGESTDT_pRank_PDEG_", PDEG, 
                "_PG1_", PA, "_tries_", n, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)



  #}  

}

#crossValidKadota(PDEG = 0.2, PA = 0.5, n = 10, iteration = 20)
#crossValidKadota(PDEG = 0.2, PA = 0.9, n = 10, iteration = 20)
#crossValidKadota(PDEG = 0.2, PA = 0.7, n = 10, iteration = 20)
#crossValidKadota(PDEG = 0.2, PA = 0.6, n = 10, iteration = 20)
#crossValidKadota(PDEG = 0.2, PA = 0.8, n = 10, iteration = 20)
#crossValidKadota(PDEG = 0.2, PA = 1,   n = 10, iteration = 20)



summaryCrossValidKadota <- function() {
  TRY <- 40
  ITE <- 20
  #PA <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
  PA <- c(0.5, 0.7, 0.9)
  REP <- c("multi", "multi", "multi", "multi", "mono", "mono", "mono")
  TAG <- c("DEGESedgeR", "DEGESTbT", "DEGESTDT", "DEGESDESeq", "DEGESDESeq", "DEGESTbT", "DEGESTDT")


  auc <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  names(auc) <- paste(REP, TAG, sep = "-")

  # read raw data.
  for (pa in 1:length(PA)) {
    for (t in 1:TRY) {
      for (tag in 1:length(TAG)) {
        if (is.null(auc[[paste(REP[tag], TAG[tag], sep = "-")]]))
          auc[[paste(REP[tag], TAG[tag], sep = "-")]] <- array(0, dim = c(length(PA), TRY, ITE + 1))
        dt <- read.table(paste(result.path, "/simResult_", REP[tag], "_Kdt_", TAG[tag],
                         "_pRank_PDEG_0.2_PG1_", PA[pa], "_tries_", t, ".txt", sep = ""), header = T)
        for (i in 1:(ncol(dt) - 2)) {
          auc[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t, i] <- AUC(rocdemo.sca(
                                truth = as.numeric(dt[, ncol(dt)] != 0), data = -dt[, i]))
        }
      }
    }
  }


  ft <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  names(ft) <- paste(REP, TAG, sep = "-")
  vl <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  names(vl) <- paste(REP, TAG, sep = "-")
  

  for (pa in 1:length(PA)) {
    for (t in 1:TRY) {
      for (tag in 1:length(TAG)) {
        if (is.null(ft[[paste(REP[tag], TAG[tag], sep = "-")]])) {
          ft[[paste(REP[tag], TAG[tag], sep = "-")]] <- array("", dim = c(length(PA), TRY))
          vl[[paste(REP[tag], TAG[tag], sep = "-")]] <- array(0, dim = c(length(PA), TRY))
          rownames(ft[[paste(REP[tag], TAG[tag], sep = "-")]]) <- PA
          rownames(vl[[paste(REP[tag], TAG[tag], sep = "-")]]) <- PA
        }
        org.auc <- auc[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t, ]
        org.auc <- org.auc[-21]
        unq.auc <- unique(org.auc)
        if (length(org.auc) == length(unq.auc)) {
          ft[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- "R"
          vl[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- NA
        } else {
          if (org.auc[length(org.auc)] == org.auc[length(org.auc) - 1]) {
            ft[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- "C"
            vl[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- 
              sum((org.auc == org.auc[length(org.auc)]) == FALSE)
          } else  {
            ft[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- "P"
            tmp.tag <- (org.auc == org.auc[length(org.auc)])
            tmp.pos <- 1:length(tmp.tag)
            tmp.pos <- tmp.pos[tmp.tag]
            vl[[paste(REP[tag], TAG[tag], sep = "-")]][pa, t] <- 
              length(tmp.tag) - max(tmp.pos[1:(length(tmp.pos) - 1)])
          }
		}
        #browser()
      }
    }
  }

  tmpl <- matrix(0, ncol = length(PA), nrow = length(ft))
  rownames(tmpl) <- names(ft)
  colnames(tmpl) <- PA
  P.num <- tmpl
  C.num <- tmpl
  R.num <- tmpl
  P.ave <- tmpl
  C.ave <- tmpl

  for (i in 1:length(ft)) {
    for (rw in 1:nrow(ft[[i]])) {
      P.tag <- (ft[[i]][rw, ] == "P")
      C.tag <- (ft[[i]][rw, ] == "C")
      R.tag <- (ft[[i]][rw, ] == "R")
      P.num[i, rw] <- sum(P.tag)
      C.num[i, rw] <- sum(C.tag)
      R.num[i, rw] <- sum(R.tag)
      P.ave[i, rw] <- mean(vl[[i]][rw, P.tag])
      C.ave[i, rw] <- mean(vl[[i]][rw, C.tag])
    }
  }
  
  table(vl[[1]][ft[[1]] == "P"])
  table(vl[[1]][ft[[1]] == "C"])
  table(vl[[2]][ft[[2]] == "P"])
  table(vl[[2]][ft[[2]] == "C"])
  table(vl[[3]][ft[[3]] == "P"])
  table(vl[[3]][ft[[3]] == "C"])
  table(vl[[4]][ft[[4]] == "P"])
  table(vl[[4]][ft[[4]] == "C"])
  table(vl[[5]][ft[[5]] == "P"])
  table(vl[[5]][ft[[5]] == "C"])
  table(vl[[6]][ft[[6]] == "P"])
  table(vl[[6]][ft[[6]] == "C"])
  table(vl[[7]][ft[[7]] == "P"])
  table(vl[[7]][ft[[7]] == "C"])
  

}









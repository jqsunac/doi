#
# This script is for analyzing tasRNA data set which is referenced from baySeq authors.
# There two arguments should set before execution.
load("./refData/rdr6_wt.RData")  # 'countData' object defained on baySeq.
result.path <- "./result_revise"        # The path for saving the results.
#
#
# Reuqirement:
#    R 2.15.2 (2012-10-26)
#    TCC 1.0.0
#    ROC 1.34.0
#    DESeq 1.10.1
#    baySeq 1.12.0
#    edgeR 3.0.4
#
#
#
# Load R libraries.
library(TCC)

# Line type, Color options.
cols <- list(TMM = "blue", TbT = "red", DES = "black")





analyzeTasRNA <- function(unique = FALSE, CPM = FALSE, floorPDEG = NULL, multi = FALSE, mono = FALSE, cycles = 10) {
  # Prepare arguments (Real data). 
  count <- cbind(as.matrix(cD@data), as.numeric(cD@annotation$TP))
  colnames(count) <- c("WT1", "WT2", "KO1", "KO2", "TP")
  tag <- "original"
  if (unique) {
    count <- unique(count)
    tag <- "unique"
  }
  if (CPM) {
    count <- sweep(count, 2, 1000000 / colSums(count), "*")
    count <- round(count)
    tag <- "rpkm"
  }
  if (!is.null(floorPDEG)) {
    tag <- "floorPDEG"
  }
  trueDEG <- count[, 5]
  count <- count[, 1:4]
  # WT1, WT2 vs KO1, KO2
  if (multi) {
  auc <- matrix(0, ncol = 2, nrow = cycles + 1)
  rownames(auc) <- paste("DEGES", c(0:cycles), sep = "")
  colnames(auc) <- c("iDEGES/edgeR", "iDEGES/TbT")
  des.fdr <- matrix(0, ncol = cycles, nrow = nrow(count))
  colnames(des.fdr) <- paste("DEGES", c(1:cycles), sep = "")
  des.potentialDEG <- matrix(0, ncol = cycles, nrow = nrow(count))
  colnames(des.potentialDEG) <- paste("DESGES", c(1:cycles), sep = "")
  tbt.fdr <- matrix(0, ncol = cycles, nrow = nrow(count))
  colnames(tbt.fdr) <- paste("TbT", c(1:cycles), sep = "")
  tbt.potentialDEG <- matrix(0, ncol = cycles, nrow = nrow(count))
  colnames(tbt.potentialDEG) <- paste("TbT", c(1:cycles), sep = "")
  norm <- c("tmm", "tmm")
  test <- c("edger", "bayseq")
  for (n in 1:length(norm)) {
    for (i in 0:cycles) {
      tcc <- TCC(count = count, group = c(2, 2))
      tcc$simulation$trueDEG <- trueDEG
      tcc <- calcNormFactors(tcc, norm.method = norm[n], test.method = test[n],
                             iteration = i, samplesize = 10000, floorPDEG = floorPDEG)
      if (n == 1 && i > 0) {
        des.fdr[, i] <- tcc$private$stat$q.value
        des.potentialDEG[, i] <- as.numeric(tcc$private$DEGES.potentialDEG != 0)
      }
      if (n == 2 && i > 0) {
        tbt.fdr[, i] <- tcc$private$stat$q.value
        tbt.potentialDEG[, i] <- as.numeric(tcc$private$DEGES.potentialDEG != 0)
      }
      tcc <- estimateDE(tcc, test.method = "edger")
      auc[i + 1, n] <- calcAUCValue(tcc)
      write.table(auc, file = paste(result.path, "/anResult_tasRNA_", tag, "_multi_AUC.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
      write.table(des.fdr, file = paste(result.path, "/anResult_tasRNA_", tag, "_multi_DEGES.FDR.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
      write.table(des.potentialDEG, file = paste(result.path, "/anResult_tasRNA_", tag, "_multi_DEGES.potentialDEG.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
      write.table(tbt.fdr, file = paste(result.path, "/anResult_tasRNA_", tag, "_multi_TbT.FDR.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
      write.table(tbt.potentialDEG, file = paste(result.path, "/anResult_tasRNA_", tag, "_multi_TbT.potentialDEG.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    }
  }
  }
  # WT1 vs KO1 and WT2 vs KO2
  if (mono) {
  norm <- c("deseq", "tmm", "tmm")
  test <- c("deseq", "deseq", "bayseq")
  auc <- matrix(0, ncol = 3 + 3, nrow = cycles + 1)
  colnames(auc) <- c("1_DEGES/deseqdeseq", "1_DEGES/tmmdeseq", "1_DEGES/tmmbayseq",
                     "2_DEGES/deseqdeseq", "2_DEGES/tmmdeseq", "2_DEGES/tmmbayseq")
  rownames(auc) <- paste("DEGES", c(0:cycles), sep = "")
  for (d in 1:2) {
    for (n in 1:length(norm)) {
      for (i in 0:cycles) {
        tryCatch({
            c <- count[, c(d, d + 2)]
            tcc <- TCC(count = c, group = c(1, 1))
            tcc$simulation$trueDEG <- trueDEG
            tcc <- calcNormFactors(tcc, norm.method = norm[n],
                                   test.method = test[n], iteration = i,
                                   floorPDEG = 0.05, samplesize = 10000)
            tcc <- estimateDE(tcc, test.method = "deseq")
            auc[i + 1, n + (d - 1) * 3] <- calcAUCValue(tcc)
          },
          error = function(e) {
            message(e)
          }, 
          finnaly = {
            auc[i + 1, n + (d - 1) * 3] <- -1
          },
          silent = TRUE
        )
        write.table(auc, file = paste(result.path, "/anResult_tasRNA_", tag, "_mono_AUC.txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
      }
    }
  }
  }
}


#analyzeTasRNA(cycles = 20, multi = TRUE)
#analyzeTasRNA(cycles = 5, mono = TRUE)



searchCyclicReasonPvalue <- function(cycles = 10) {
  count <- cbind(cD@data, cD@annotation$TP)
  colnames(count) <- c("WT1", "WT2", "KO1", "KO2", "TP")
  tp <- as.numeric(cD@annotation$TP)
  a.count <- count[, 1:4]
  a.tp <- count[, 5]
  unique.count <- unique(count)
  b.count <- unique.count[, 1:4]
  b.tp <- unique.count[, 5]
  c.count <- rbind(b.count, matrix(rep(c(1,0,0,0), times = 10000), ncol = 4, byrow = TRUE))
  c.tp <- c(b.tp, rep(0, length = 10000))
  
  norm <- c("tmm", "tmm", "tmm")
  test <- c("edger", "bayseq", "deseq")

  auc <- matrix(0, ncol = 6, nrow = cycles + 1)
  colnames(auc) <- c("tmm-edger-a", "tmm-edger-c", "tmm-bayseq-a", "tmm-bayseq-c",
                     "tmm-deseq-a", "tmm-deseq-c")

  tbt.lock <- FALSE

  for (i in 0:cycles) {
    for (m in 1:length(norm)) {
      if (m == 2 && tbt.lock)
        next
      tcc <- TCC(count = a.count, group = c(2, 2))
      tcc$simulation$trueDEG <- as.numeric(a.tp)
      tcc <- calcNormFactors(tcc, norm.method = norm[m], test.method = test[m],
                             iteration = i, samplesize = 10000)
      tcc <- estimateDE(tcc, "edger")
      auc[i + 1, 1 + 2 * (m - 1)] <- calcAUCValue(tcc)
      write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_pval.txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      if (i >= 1 && auc[i, 3] == auc[i + 1, 3])
        tbt.lock <- TRUE
    }
  }

  tbt.lock <- FALSE
  for (i in 0:cycles) {
    for (m in 1:length(norm)) {
      if (m == 2 && tbt.lock)
        next
      tcc <- TCC(count = c.count, group = c(2, 2))
      tcc$simulation$trueDEG <- as.numeric(c.tp)
      tcc <- calcNormFactors(tcc, norm.method = norm[m], test.method = test[m],
                             iteration = i, samplesize = 10000)
      tcc <- estimateDE(tcc, "edger")
      auc[i + 1, 2 * m] <- calcAUCValue(tcc)
      write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_pval.txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      if (i >= 1 && auc[i, 4] == auc[i + 1, 4])
        tbt.lock <- TRUE
    }
  }
}

#searchCyclicReasonPvalue(cycles = 50)


searchCyclicReason <- function(cycles = 10, dataset = "a") {
  count <- cbind(cD@data, cD@annotation$TP)
  colnames(count) <- c("WT1", "WT2", "KO1", "KO2", "TP")
  tp <- as.numeric(cD@annotation$TP)

  a.count <- count[, 1:4]
  a.tp <- count[, 5]

  unique.count <- unique(count)
  b.count <- unique.count[, 1:4]
  b.tp <- unique.count[, 5]

  c.count <- rbind(b.count, matrix(rep(c(1,0,0,0), times = 10000), ncol = 4, byrow = TRUE))
  c.tp <- c(b.tp, rep(0, length = 10000))
  d.count <- rbind(c.count, matrix(rep(c(0,1,0,0), times = 20000), ncol = 4, byrow = TRUE))
  d.tp <- c(c.tp, rep(0, length = 20000))
  e.count <- rbind(d.count, matrix(rep(c(0,0,0,1), times = 10000), ncol = 4, byrow = TRUE))
  e.tp <- c(d.tp, rep(0, length = 10000))

  s.count <- rbind(b.count, matrix(rep(c(32,22,30,25), times = 10000), ncol = 4, byrow = TRUE))
  s.tp <- c(b.tp, rep(0, length = 10000))
  t.count <- rbind(s.count, matrix(rep(c(58,49,52,55), times = 20000), ncol = 4, byrow = TRUE))
  t.tp <- c(s.tp, rep(0, length = 20000))


  auc <- matrix(0, ncol = 9, nrow = cycles + 1)
  rownames(auc) <- paste("DEGES", c(0:cycles), sep = "")
  colnames(auc) <- c("a:original", "b:unique", "c", "d", "e", "s", "t", "a:DESeq", "c:DESeq")


  if (dataset == "a") {
  for (i in 0:cycles) {
    tcc <- TCC(count = a.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(a.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 1] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
  if (dataset == "b") {
  for (i in 0:cycles) {
    tcc <- TCC(count = b.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(b.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 2] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
  if (dataset == "c") {
  for (i in 0:cycles) {
    tcc <- TCC(count = c.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(c.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 3] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
  if (dataset == "d") {
  for (i in 0:cycles) {
    tcc <- TCC(count = d.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(d.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 4] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
 if (dataset == "e") {
  for (i in 0:cycles) {
    tcc <- TCC(count = e.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(e.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 5] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
  if (dataset == "s") {
  for (i in 0:cycles) {
    tcc <- TCC(count = s.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(s.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 6] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
  if (dataset == "t") {
  for (i in 0:cycles) {
    tcc <- TCC(count = t.count, group = c(2, 2))
    tcc$simulation$trueDEG <- as.numeric(t.tp)
    tcc <- calcNormFactors(tcc, "tmm", "edger", iteration = i)
    tcc <- estimateDE(tcc, "edger")
    auc[i + 1, 7] <- calcAUCValue(tcc)
    write.table(auc, file = paste(result.path, "/anResult_tasRNA_dataset_", dataset, ".txt", sep = ""),
          sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  }
}

unionCyclicReason <- function() {
  a <- read.table(paste(result.path, "/anResult_tasRNA_dataset_a.txt", sep = ""), header = TRUE)
  b <- read.table(paste(result.path, "/anResult_tasRNA_dataset_b.txt", sep = ""), header = TRUE)
  c <- read.table(paste(result.path, "/anResult_tasRNA_dataset_c.txt", sep = ""), header = TRUE)
  d <- read.table(paste(result.path, "/anResult_tasRNA_dataset_d.txt", sep = ""), header = TRUE)
  e <- read.table(paste(result.path, "/anResult_tasRNA_dataset_e.txt", sep = ""), header = TRUE)
  s <- read.table(paste(result.path, "/anResult_tasRNA_dataset_s.txt", sep = ""), header = TRUE)
  t <- read.table(paste(result.path, "/anResult_tasRNA_dataset_t.txt", sep = ""), header = TRUE)
  df <- data.frame(
    a = a[, 1], b = b[, 2], c = c[, 3], d = d[, 4],
    e = e[, 5], s = s[, 6], t = t[, 7]
  )
  write.table(df, file = paste(result.path, "/anResult_tasRNA_dataset_all.txt", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#searchCyclicReason(cycles = 100, dataset = "a")
#searchCyclicReason(cycles = 100, dataset = "b")
#searchCyclicReason(cycles = 100, dataset = "c")
#searchCyclicReason(cycles = 100, dataset = "d")
#searchCyclicReason(cycles = 100, dataset = "e")
#searchCyclicReason(cycles = 100, dataset = "s")
#searchCyclicReason(cycles = 100, dataset = "t")
#unionCyclicReason()


getTable6 <- function() {
  potentialDEG <- read.table(paste(result.path,
                             "/anResult_tasRNA_original_multi_DEGES.potentialDEG.txt",
                             sep = ""), header = TRUE)
  count <- cbind(cD@data)
  colnames(count) <- c("WT1", "WT2", "KO1", "KO2")
  cy <- 9:15
  x <- potentialDEG[, cy]

  sum(rowSums(x) == 7)
  sum(rowSums(x) > 0)
  sum(rowSums(x) == 6 & x[, 1] != 1)
  sum(rowSums(x) == 6 & x[, 2] != 1)
  sum(rowSums(x) == 6 & x[, 3] != 1)
  sum(rowSums(x) == 6 & x[, 4] != 1)
  sum(rowSums(x) == 6 & x[, 5] != 1)
  sum(rowSums(x) == 6 & x[, 6] != 1)
  sum(rowSums(x) == 6 & x[, 7] != 1)
  sum(rowSums(x) == 1 & x[, 1] == 1)
  sum(rowSums(x) == 1 & x[, 2] == 1)
  sum(rowSums(x) == 1 & x[, 3] == 1)
  sum(rowSums(x) == 1 & x[, 4] == 1)
  sum(rowSums(x) == 1 & x[, 5] == 1)
  sum(rowSums(x) == 1 & x[, 6] == 1)
  sum(rowSums(x) == 1 & x[, 7] == 1)

  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 2] != 1)
  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 3] != 1)
  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 4] != 1)
  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 5] != 1)
  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 6] != 1)
  sum(rowSums(x) == 5 & x[, 1] != 1 & x[, 7] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 1] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 2] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 3] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 4] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 5] != 1)
  sum(rowSums(x) == 5 & x[, 6] != 1 & x[, 7] != 1)

  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 2] == 1)
  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 3] == 1)
  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 4] == 1)
  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 5] == 1)
  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 6] == 1)
  sum(rowSums(x) == 2 & x[, 1] == 1 & x[, 7] == 1)
  sum(rowSums(x) == 2 & x[, 2] == 1 & x[, 3] == 1)
  sum(rowSums(x) == 2 & x[, 2] == 1 & x[, 4] == 1)
  sum(rowSums(x) == 2 & x[, 2] == 1 & x[, 5] == 1)
  sum(rowSums(x) == 2 & x[, 2] == 1 & x[, 6] == 1)
  sum(rowSums(x) == 2 & x[, 2] == 1 & x[, 7] == 1)
  sum(rowSums(x) == 2 & x[, 3] == 1 & x[, 4] == 1)
  sum(rowSums(x) == 2 & x[, 3] == 1 & x[, 5] == 1)
  sum(rowSums(x) == 2 & x[, 3] == 1 & x[, 6] == 1)
  sum(rowSums(x) == 2 & x[, 3] == 1 & x[, 7] == 1)
  sum(rowSums(x) == 2 & x[, 4] == 1 & x[, 5] == 1)
  sum(rowSums(x) == 2 & x[, 4] == 1 & x[, 6] == 1)
  sum(rowSums(x) == 2 & x[, 4] == 1 & x[, 7] == 1)
  sum(rowSums(x) == 2 & x[, 5] == 1 & x[, 6] == 1)
  sum(rowSums(x) == 2 & x[, 5] == 1 & x[, 7] == 1)
  sum(rowSums(x) == 2 & x[, 6] == 1 & x[, 7] == 1)


  count.ptPDEG <- count[(rowSums(x) > 0), ]

  z <- count[(rowSums(x) == 7), ]
  w <- count[(rowSums(x) > 0) & (rowSums(x) != 7), ]


  nf <- matrix(0, ncol = length(cy), nrow = ncol(count))
  for (i in cy) {
    tcc <- new("TCC", count, group = c(1, 1, 2, 2))
    tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = i)
    nf[, i - 8] <- tcc$norm.factors
  }

  t <- combn(ncol(x), 2, function(ij) {
    i <- ij[1]
    j <- ij[2]
	#cor
    #return(cor(x[, i], x[, j]))
    #eculid
	dis <- sqrt(sum( (x[, i] - x[, j]) ^ 2 ))
	#browser()
    return(dis)
	})
  cb <- combn(ncol(x), 2)
  d <- matrix(0, ncol = ncol(x), nrow = ncol(x))
  colnames(d) <- cy
  colnames(d) <- cy
  for (col in 1:ncol(cb)) {
    d[cb[1, col], cb[2, col]] <- 1 - t[col]
  }
  d2 <- dist(t(d))
  h <- hclust(d2)
  png(paste(result.path, "/smallRNA_clust.png", sep = ""), 400, 500)
  plot(h, hang = -1)
  dev.off()

  
  
  
  
}






plotFigure2 <- function(width = 720, height = 330, unique = FALSE, CPM = FALSE, floorPDEG = FALSE) {
  tag <- "original"
  if (unique)
    tag <- "unique"
  if (CPM)
    tag <- "rpkm"
  if(floorPDEG)
    tag <- "floorPDEG"

  # Figure 2-a
  png(paste(result.path, "/Figure-2a-", tag, ".png", sep = ""), width = width, height = height)
  auc <- read.table(paste(result.path, "/anResult_tasRNA_", tag, "_multi_AUC.txt", sep = ""), header = TRUE)
  auc <- auc * 100
  if (tag == "original")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(60, 70), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "unique")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(92, 93), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "floorPDEG")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(60, 75), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "rpkm")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(70, 80), xlab = "iteration", ylab = "AUC (%)")
  grid()
  points(0:10, auc[1:11, 1], col = cols$DES)   # iDES/TMM
  lines(0:10, auc[1:11, 1], col = cols$DES)    # iDEs/TMM
  points(0:10, auc[1:11, 2], col = cols$TbT)   # iDES/TbT
  lines(0:10, auc[1:11, 2], col = cols$TbT)    # iDEs/TbT
  legend("bottomright", legend = c("iDEGES/TbT", "iDEGES/edgeR"), col = c(cols$TbT, cols$DES),
         lty = 1, pch = 1)
  dev.off()

  # Figure 2-c
  png(paste(result.path, "/Figure-2c-", tag, ".png", sep = ""), width = width, height = height)
  auc <- read.table(paste(result.path, "/anResult_tasRNA_", tag, "_multi_AUC.txt", sep = ""), header = TRUE)
  auc <- auc * 100
  if (tag == "original")
    plot(0, 0, type = "n", xlim = c(0, 50), ylim = c(60, 70), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "unique")
    plot(0, 0, type = "n", xlim = c(0, 50), ylim = c(92, 93), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "floorPDEG")
    plot(0, 0, type = "n", xlim = c(0, 50), ylim = c(60, 75), xlab = "iteration", ylab = "AUC (%)")
  if (tag == "rpkm")
    plot(0, 0, type = "n", xlim = c(0, 50), ylim = c(70, 80), xlab = "iteration", ylab = "AUC (%)")
  grid()
  if (tag == "original") {
    points(0:(nrow(auc) - 1), auc[, 1], col = cols$DES)   # iDES/TMM
    lines(0:(nrow(auc) - 1), auc[, 1], col = cols$DES)    # iDEs/TMM
    points(0:10, auc[1:11, 2], col = cols$TbT)   # iDES/TbT
    lines(0:10, auc[1:11, 2], col = cols$TbT)    # iDEs/TbT
  } else {
    points(0:(nrow(auc) - 1), auc[, 1], col = cols$DES)   # iDES/TMM
    lines(0:(nrow(auc) - 1), auc[, 1], col = cols$DES)    # iDEs/TMM
    points(0:(nrow(auc) - 1), auc[, 2], col = cols$TbT)   # iDES/TbT
    lines(0:(nrow(auc) - 1), auc[, 2], col = cols$TbT)    # iDEs/TbT
	}
  legend("bottomright", legend = c("iDEGES/TbT", "iDEGES/edgeR"), col = c(cols$TbT, cols$DES),
         lty = 1, pch = 1)
  dev.off()

  # Figure 2-b
  png(paste(result.path, "/Figure-2b-", tag, ".png", sep = ""), width = width, height = height)
  des.fdr <- read.table(paste(result.path, "/anResult_tasRNA_", tag, "_multi_DEGES.FDR.txt", sep = ""), header = TRUE)
  des.pDEG <- read.table(paste(result.path, "/anResult_tasRNA_", tag, "_multi_DEGES.potentialDEG.txt", sep = ""), header = TRUE)
  tbt.pDEG <- read.table(paste(result.path, "/anResult_tasRNA_", tag, "_multi_TbT.potentialDEG.txt", sep = ""), header = TRUE)
  des.fdr <- des.fdr[, 1:10] 
  des.pDEG <- des.pDEG[, 1:10] 
  tbt.pDEG <- tbt.pDEG[, 1:10]
  labels <- c("iDEGES/TbT", "iDEGES/edgeR (used)", "iDEGES/edgeR (FDR)")
  if (tag == "original")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 7), xlab = "iteration", ylab = "PDEG (%)")
  if (tag == "unique" || tag == "rpkm")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), xlab = "iteration", ylab = "PDEG (%)")
  if (tag == "floorPDEG")
    plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 12), xlab = "iteration", ylab = "PDEG (%)")
  grid()
  points(1:10, colSums(des.fdr < 0.1) * 100 / nrow(des.fdr), col = cols$DES, pch = 1) # iDES/TMM
  lines(1:10, colSums(des.fdr < 0.1) * 100 / nrow(des.fdr), col = cols$DES, lty = 2)  # iDES/TMM
  points(1:10, colSums(des.pDEG) * 100 / nrow(des.pDEG), col = cols$DES, pch = 1)  # iDES/TMM (used)
  lines(1:10, colSums(des.pDEG) * 100 / nrow(des.pDEG), col = cols$DES, lty = 1)   # iDES/TMM (used)
  points(1:10, colSums(tbt.pDEG) * 100 / nrow(tbt.pDEG), col = cols$TbT, pch = 1)  # iDES/TbT (used)
  lines(1:10, colSums(tbt.pDEG) * 100 / nrow(tbt.pDEG), col = cols$TbT, lty = 1)   # iDES/TbT (used)
  legend("right", legend = labels, pch = 1, lty = c(1, 1, 2), col = c(cols$TbT, cols$DES, cols$DES))
  dev.off()
}



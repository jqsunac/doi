library(TCC)
result.path <- "./result_revise"


mainSimCore <- function(PDEG = NULL, PA = NULL, SEED = NULL, rep = NULL) {
  type.tag <- "ALL"
  samplesize <- 10000
  norm.norm <- NULL
  norm.test <- NULL
  norm.iter <- NULL
  test.test <- NULL
  if (rep == "multi") {
    norm.norm <- c("tmm",  "tmm",    "tmm",   "tmm",   "deseq", "tmm",   "deseq",  "tmm")
    norm.test <- c("none", "bayseq", "edger", "deseq", "deseq", "edger", "deseq", "deseq")
    norm.iter <- c(0,      1,        1,       1,       1,       3,       3,       3)
    test.test <- c("edger", "deseq", "bayseq")
  } 
  if (rep == "mono") {
    norm.norm <- c("deseq", "deseq", "tmm",   "tmm",    "deseq", "tmm")
    norm.test <- c("none",  "deseq", "deseq", "bayseq", "deseq", "deseq")
    norm.iter <- c(0,       1,       1,       1,         3,      3)
    test.test <- c("deseq", "bayseq")
  }

  set.seed(SEED)
  if (rep == "multi") {
    tcc <- simulateReadCounts(Ngene = 10000, PDEG = PDEG, DEG.assign = c(PA, 1 - PA))
  }
  if (rep == "mono") {
    tcc <- simulateReadCounts(Ngene = 10000, PDEG = PDEG, DEG.assign = c(PA, 1 - PA),
                              replicates = c(1, 1))
  }
  for (n in 1:length(norm.norm)) {
    # normalization
    tcc.n <- calcNormFactors(tcc, norm.method = norm.norm[n], test.method = norm.test[n],
                             iteration = norm.iter[n], samplesize = samplesize)
    # write execution time
    write.table(as.numeric(tcc.n$DEGES$execution.time), paste(result.path, "/sim_SIM", rep, "TIME", 
                norm.norm[n], norm.test[n], norm.iter[n], "none", 
                "PDEG", PDEG, "PA", PA, "SEED", SEED, ".txt", sep = "_"),
                col.names = TRUE, row.names = FALSE)
    for (t in 1:length(test.test)) {
      # identifying
      tcc.t <- estimateDE(tcc.n, test.method = test.test[t], samplesize = samplesize)
      pdeg <- rep(0, length = nrow(tcc.t$count))
      pval <- rep(0, length = nrow(tcc.t$count))
      if (length(tcc.t$private$debug$potentialDEG) > 0) {
        pdeg <- tcc.t$private$debug$potentialDEG[[length(tcc.t$private$debug$potentialDEG)]]
        pval <- tcc.t$private$debug$potentialPval[[length(tcc.t$private$debug$potentialPval)]]
      }
      df <- data.frame(
        TRUETH = as.numeric(tcc$simulation$trueDEG != 0),
        ESTIMATED = as.numeric(tcc.t$estimatedDEG != 0),
        PDEG = pdeg,
        PVAL = pval
      )
      # write statistics
      write.table(df, paste(result.path, "/sim_SIM", rep, type.tag, 
                  norm.norm[n], norm.test[n], norm.iter[n], test.test[t], 
                  "PDEG", PDEG, "PA", PA, "SEED", SEED, ".txt", sep = "_"),
                  col.names = TRUE, row.names = FALSE)
    }
  }
}


simMultiple <- function() {
  PDEG <- c(0.05, 0.15, 0.25)
  PA <- c(0.50, 0.70, 0.90)
  SEED <- 1:100
  for (a in 1:length(PDEG)) {
    for (d in 1:length(PA)) {
      for (s in 1:length(SEED)) {
        mainSimCore(PDEG = PDEG[d], PA = PA[a], SEED = SEED[s], rep = "multi")
      }
    }
  }
}

simMono <- function(SEED) {
  PDEG <- c(0.05, 0.15, 0.25)
  PA <- c(0.50, 0.70, 0.90)
  for (a in 1:length(PDEG)) {
    for (d in 1:length(PA)) {
      for (s in 1:length(SEED)) {
        mainSimCore(PDEG = PDEG[d], PA = PA[a], SEED = SEED[s], rep = "mono")
      }
    }
  }
}

simMono(SEED = 20)



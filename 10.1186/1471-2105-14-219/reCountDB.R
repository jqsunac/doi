#
data.path <- "./refData"
result.path <- "./result_revise"        # The path for saving the results.
#
#
# Load R libraries.
library(TCC)



loopDES <- function(ite = -1, count, group, filecode = "",
                      unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- count[(rowSums(count) > 0), ]
  tag <- "original"
  if (unique) {
    count <- unique(count)
    tag <- "unique"
  }
  if (RPM) {
    count <- round(sweep(count, 2, 1000000 / colSums(count), "*"))
    tag <- "rpm"
  }

  if (ncol(count) == 2) {
    norm <- c("deseq", "tmm", "tmm")
    test <- c("deseq",  "deseq", "bayseq")
  } else {
    norm <- c("tmm", "tmm", "deseq", "tmm")
    test <- c("edger",  "deseq", "deseq", "bayseq")
  }  

  for (n in 1:length(norm)) {
    tcc <- new("TCC", count, group)
    tcc <- calcNormFactors(tcc, norm.method = norm[n], test.method = test[n],
                           iteration = ite, samplesize = 10000, floorPDEG = floorPDEG)
    pRank <- matrix(0, ncol = ite, nrow = nrow(count))
    pType <- rep("NONE", length = ite)
    pInput <- rep(0, length = ite)
    pPDEG <- rep(0, length = ite)
    for (i in 1:ite) {
      pRank[, i] <- tcc$private$debug$potentialPval[[i]]
      pType[i] <- tcc$private$debug$potentialThreshold[[i]]$type
      pInput[i] <- tcc$private$debug$potentialThreshold[[i]]$input
      pPDEG[i] <- tcc$private$debug$potentialThreshold[[i]]$PDEG
    }
    pThre <- cbind(type = pType, input = pInput, PDEG = pPDEG)
    if (is.null(floorPDEG)) {
      flp <- ""
    } else {
      flp <- floorPDEG
    }
    write.table(pRank, file = paste(result.path, "/anResult_recount_", filecode,
                "_", norm[n], "_", test[n],
                "_pRank_", tag, "_floorPDEG_", flp, 
                "_.txt", sep = ""), col.names = FALSE, row.names = FALSE)
    write.table(pThre, file = paste(result.path, "/anResult_recount_", filecode,
                "_", norm[n], "_", test[n],
                "_pThre_", tag, "_floorPDEG_", flp, 
                "_.txt", sep = ""), col.names = FALSE, row.names = FALSE)
  }
}




analyzeGilad <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/gilad_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1,1,1,2,2,2)
  loopDES(n, count, group, "gilad", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}

analyzeMaqc <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/maqc_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1,1,1,1,1,1,1,2,2,2,2,2,2,2)
  loopDES(n, count, group, "maqc", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}
analyzeMaqcP <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/maqcpooled_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1,2)
  loopDES(n, count, group, "maqcpooled", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}


analyzeMontpick <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/montpick_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,
             2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2)
  loopDES(n, count, group, "montpick", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}


analyzeKatzmouse  <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/katzmouse_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  count <- count[, c(1, 3, 2, 4)]
  group <- c(1, 1, 2, 2)
  loopDES(n, count, group, "katzmouse", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}


analyzeBottomly  <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/bottomly_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2)
  loopDES(n, count, group, "bottomly", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}

analyzeSultan <- function(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL) {
  count <- read.table(paste(data.path, "/sultan_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  group <- c(1, 1, 2, 2)
  loopDES(n, count, group, "sultan", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}

analyzeRdr6 <- function(n, unique = FALSE, RPM = RPM, floorPDEG = floorPDEG) {
  load(paste(data.path, "/rdr6_wt.RData", sep = ""))
  count <- as.matrix(cD@data)
  group <- c(1, 1, 2, 2)
  loopDES(n, count, group, "rdr6", unique = unique, RPM = RPM, floorPDEG = floorPDEG)
}



getReplicates <- function() {
  recount.files <- c(
    paste(data.path, "/bottomly_count_table.txt", sep = ""),
    paste(data.path, "/gilad_count_table.txt", sep = ""),
    paste(data.path, "/katzmouse_count_table.txt", sep = ""),
    paste(data.path, "/maqc_count_table.txt", sep = ""),
    paste(data.path, "/sultan_count_table.txt", sep = "")
  )

  for (i in 1:length(recount.files)) {
    count <- read.table(recount.files[i], header = TRUE)
    count <- count[, - 1]
    count <- count[(sum(rowSums(count)) > 0), ]
    cat(recount.files[i])
    cat("\n")
    cat("original  ")
    cat(dim(count))
    cat("   ")
    cat(sum(colSums(count)) / nrow(count))
    unique <- unique(count)
    cat("\n")
    cat("unique  ")
    cat(dim(unique))
    cat("   ")
    cat(sum(colSums(unique)) / nrow(unique))
    cat("\n")
  }

}


n <- 30

#analyzeGilad(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL)
#analyzeMaqc(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL)
#analyzeMaqcP(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL)
#analyzeKatzmouse(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL)
#analyzeSultan(n, unique = FALSE, RPM = FALSE, floorPDEG = NULL)

#analyzeGilad(n, unique = TRUE, RPM = FALSE, floorPDEG = NULL)
#analyzeMaqc(n, unique = TRUE, RPM = FALSE, floorPDEG = NULL)
#analyzeMaqcP(n, unique = TRUE, RPM = FALSE, floorPDEG = NULL)
#analyzeKatzmouse(n, unique = TRUE, RPM = FALSE, floorPDEG = NULL)
#analyzeSultan(n, unique = TRUE, RPM = FALSE, floorPDEG = NULL)


#analyzeGilad(n, unique = FALSE, RPM = FALSE, floorPDEG = 0.1)
#analyzeMaqc(n, unique = FALSE, RPM = FALSE, floorPDEG = 0.1)
#analyzeMaqcP(n, unique = FALSE, RPM = FALSE, floorPDEG = 0.1)
#analyzeKatzmouse(n, unique = FALSE, RPM = FALSE, floorPDEG =0.1)
#analyzeSultan(n, unique = FALSE, RPM = FALSE, floorPDEG = 0.1)

#analyzeGilad(n, unique = TRUE, RPM = FALSE, floorPDEG = 0.1)
#analyzeMaqc(n, unique = TRUE, RPM = FALSE, floorPDEG = 0.1)
#analyzeMaqcP(n, unique = TRUE, RPM = FALSE, floorPDEG = 0.1)
#analyzeKatzmouse(n, unique = TRUE, RPM = FALSE, floorPDEG = 0.1)
#analyzeSultan(n, unique = TRUE, RPM = FALSE, floorPDEG = 0.1)







  



result.path <- "result_revise"
chkFeatures <- function(pDEGs) {
  ft <- list()
  ft$type <- "NONE"
  ft$val <- -1

  val <- rep(0, length = ncol(pDEGs))
  for (i in (length(val)):1) {
    vec <- pDEGs[, i] - pDEGs[, length(val)]
    val[i] <- sum(vec != 0)
  }
  zero.pos <- 0
  zero.1st <- TRUE
  for (i in 1:(length(val))) {
    if (val[i] == 0) {
      if (zero.1st) {
        zero.pos <- i
        zero.1st <- FALSE
        ft$type <- "R"
        ft$val <- 0
      } else {
        if (zero.pos + 1 == i) {
          ft$type <- "C"
          ft$val <- zero.pos
          break
        } else {
          ft$type <- "P"
          ft$val <- i - zero.pos
          break
        }
      }
    } else {
      ft$type <- "R"
      ft$val <- 0
    }
  }
  return (ft)
}

anaLoops <- function() {
  norm <- c("tmm", "tmm", "deseq", "tmm")
  test <- c("edger", "bayseq", "deseq", "deseq")
  type <- c("original", "unique")
  file.tag <- c("rdr6", "gilad", "maqc", "maqcpooled", "katzmouse", "sultan")
  df05 <- data.frame(
    featureType = rep("N", length = length(file.tag)),
    featureVal = rep(0, length = length(file.tag)),
    thresholdType = rep("N", length = length(file.tag)),
    thresholdVal = rep(0, length = length(file.tag))
  )
  rownames(df05) <- file.tag
  df10 <- data.frame(
    featureType = rep("N", length = length(file.tag)),
    featureVal = rep(0, length = length(file.tag)),
    thresholdType = rep("N", length = length(file.tag)),
    thresholdVal = rep(0, length = length(file.tag))
  )
  rownames(df10) <- file.tag
  df05 <- data.frame(
    DEGESedger_C = rep(-1, length = length(file.tag) * length(type)),
    DEGESedger_P = rep(-1, length = length(file.tag) * length(type)),
    DEGEStbt_C = rep(-1, length = length(file.tag) * length(type)),
    DEGEStbt_P = rep(-1, length = length(file.tag) * length(type)),
    DEGESdeseq_C = rep(-1, length = length(file.tag) * length(type)),
    DEGESdeseq_P = rep(-1, length = length(file.tag) * length(type)),
    DEGEStdt_C = rep(-1, length = length(file.tag) * length(type)),
    DEGEStdt_P = rep(-1, length = length(file.tag) * length(type))
  )
  df10 <- data.frame(
    DEGESedger_C = rep(-1, length = length(file.tag) * length(type)),
    DEGESedger_P = rep(-1, length = length(file.tag) * length(type)),
    DEGEStbt_C = rep(-1, length = length(file.tag) * length(type)),
    DEGEStbt_P = rep(-1, length = length(file.tag) * length(type)),
    DEGESdeseq_C = rep(-1, length = length(file.tag) * length(type)),
    DEGESdeseq_P = rep(-1, length = length(file.tag) * length(type)),
    DEGEStdt_C = rep(-1, length = length(file.tag) * length(type)),
    DEGEStdt_P = rep(-1, length = length(file.tag) * length(type))
  )

  type <- "original"
  for (f in 1:length(file.tag)) {
    for (n in 1:length(norm)) {
      for (t in 1:length(type)) {
        e <- try(
          p.deg <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                         norm[n], "_", test[n], "_pDEG_", type[t], "_floorPDEG__.txt", sep = ""), header = FALSE)
        )
        if (class(e) == "try-error")  {
          df05[f * 2 - (2 - t), n * 2 - 1] <- -2
          df05[f * 2 - (2 - t), n * 2] <- -2
          next
        }
        ft <- chkFeatures(p.deg)
        if (ft$type == "C") {
          df05[f * 2 - (2 - t), n * 2 - 1] <- ft$val
          df05[f * 2 - (2 - t), n * 2] <- 0
        }
        if (ft$type == "P") {
          df05[f * 2 - (2 - t), n * 2 - 1] <- 0
          df05[f * 2 - (2 - t), n * 2] <- ft$val
        }
      }
    }
  }
  for (f in 1:length(file.tag)) {
    for (n in 1:length(norm)) {
      for (t in 1:length(type)) {
        e <- try(
          p.deg <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                         norm[n], "_", test[n], "_pDEG_", type[t], "_floorPDEG_0.1_.txt", sep = ""), header = FALSE)
        )
        if (class(e) == "try-error")  {
          df10[f * 2 - (2 - t), n * 2 - 1] <- -2
          df10[f * 2 - (2 - t), n * 2] <- -2
          next
        }
        ft <- chkFeatures(p.deg)
        if (ft$type == "C") {
          df10[f * 2 - (2 - t), n * 2 - 1] <- ft$val
          df10[f * 2 - (2 - t), n * 2] <- 0
        }
        if (ft$type == "P") {
          df10[f * 2 - (2 - t), n * 2 - 1] <- 0
          df10[f * 2 - (2 - t), n * 2] <- ft$val
        }
      }
    }
  }

  
  prepdeg <- matrix(0, nrow = 30, ncol = length(file.tag) * length(norm))
  pdeg <- matrix(0, nrow = 30, ncol = length(file.tag) * length(norm))
  pdeg10 <- matrix(0, nrow = 30, ncol = length(file.tag) * length(norm))
  type <- c("original")
  for (f in 1:length(file.tag)) {
    for (n in 1:length(norm)) {
      for (t in 1:length(type)) {
        e <- try(
          p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                         norm[n], "_", test[n], "_pThre_", type[t], "_floorPDEG__.txt", sep = ""), header = FALSE)
        )
        if (class(e) == "try-error") {
          pdeg[, (f - 1) * length(norm) + n] <- -2
        } else {
          pdeg[, (f - 1) * length(norm) + n] <- p.rank[1:30, 3]
        }
        e <- try(
          p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                         norm[n], "_", test[n], "_pThre_", type[t], "_floorPDEG_0.1_.txt", sep = ""), header = FALSE)
        )
        if (class(e) == "try-error") {
          pdeg10[, (f - 1) * length(norm) + n] <- -2
        } else {
          pdeg10[, (f - 1) * length(norm) + n] <- p.rank[1:30, 3]
        }
        if (norm[n] == "tmm" && test[n] == "bayseq") {
          e <- try(
            p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                           norm[n], "_", test[n], "_pThre_", type[t], "_floorPDEG__.txt", sep = ""), header = FALSE)
          )
          if (class(e) == "try-error") {
            prepdeg[, (f - 1) * length(norm) + n] <- -2
          } else {
            prepdeg[, (f - 1) * length(norm) + n] <- p.rank[1:30, 2]
          }
        } else {
          e <- try(
            p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                           norm[n], "_", test[n], "_pRank_", type[t], "_floorPDEG__.txt", sep = ""), header = FALSE)
          )
          if (class(e) == "try-error") {
            prepdeg[, (f - 1) * length(norm) + n] <- -2
          } else {
            q <- apply(p.rank, 2, function(i) {
                         qval <- p.adjust(i, method = "BH")
                         return(sum(qval < 0.1) / length(qval))
                       })
            prepdeg[, (f - 1) * length(norm) + n] <- q[1:30]
          }
        }
      }
    }
  }

  write.table(pdeg, file = "recount.db.pdeg.txt", col.names = F, row.names = F)
  write.table(pdeg10, file = "recount.db.pdeg10.txt", col.names = F, row.names = F)
  write.table(prepdeg, file = "recount.db.prepdeg.txt", col.names = F, row.names = F)


  # fig
  y <- array(0, dim = c(length(file.tag), length(norm), 30))
  for (f in 1:length(file.tag)) {
    for (n in 1:length(norm)) {
      type <- c("original")
      #type <- c("unique")
      for (t in 1:length(type)) {
        e <- try(
          p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
                         norm[n], "_", test[n], "_pRank_", type[t], "_.txt", sep = ""), header = FALSE)
          #p.rank <- read.table(paste(result.path, "/anResult_recount_", file.tag[f], "_",
          #               norm[n], "_", test[n], "_pRank_", type[t], "_floorPDEG_0.1_.txt", sep = ""), header = FALSE)
        )
        if (class(e) == "try-error")  {
          y[f, n, ] <- 0
          next
        }
        for (cols in 2:30) {
          y[f, n, cols] <- cor(p.rank[, 1], p.rank[, cols])
        }
      }
    }
  }
  for (f in 1:length(file.tag)) {
    png(paste("recount_", file.tag[f], "_original_floorPDEG_0.05.png"), 350, 320)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(4, 4, 1, 1))
    plot(0, 0, type = "n", xlab = "iteration", ylab = "correlations", xlim = c(1, 30), ylim = c(0.91, 1.00))
    grid()
    cols <- c("black", "red", "lightblue", "blue")
    for (n in 1:length(norm)) {
      points(2:30, y[f, n, 2:30], col = cols[n])
      lines(2:30, y[f, n, 2:30], col = cols[n])
    }
    legend("bottomright", ncol = 2,
         legend = c("iDEGES/TbT", "iDEGES/edgeR", "iDEGES/TDT", "iDEGES/DESeq"),
         col = c("red", "black", "blue", "lightblue"), pch = 1, lty = 1, cex = 0.8)
    dev.off()
  }
}









  
giladVector <- function() {
  count <- read.table(paste(data.path, "/gilad_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  count <- count[(rowSums(count) > 0), ]
  group <- c(1, 1, 1, 2, 2, 2)

  norm.norm <- c("tmm", "tmm", "deseq")
  norm.test <- c("edger", "deseq", "deseq")
  norm.iter <- c(5, 5, 5)
  #test <- "edger"
  #test <- "bayseq"
  test <- "deseq"
  
  for (n in 1:length(norm.norm)) {
    ranks <- matrix(0, ncol = norm.iter[n] + 1, nrow = nrow(count))
    tcc <- new("TCC", count, group = group)
    tcc <- calcNormFactors(tcc, norm.method = norm.norm[n],
                           test.method = norm.test[n], iteration = 0)
    tcc <- estimateDE(tcc, test.method = test)
    ranks[, 1] <- tcc$stat$rank
    for (i in 1:norm.iter[n]) {
      tcc <- calcNormFactors(tcc, norm.method = norm.norm[n],
                           test.method = norm.test[n], increment = TRUE)
      tcc <- estimateDE(tcc, test.method = test)
      ranks[, i + 1] <- tcc$stat$rank
      write.table(ranks, paste("result_revise/gilad_RankVec", norm.norm[n], norm.test[n],
                  test, sep = "_"), row.names = FALSE, col.names = FALSE)
    }
  }
}


corLines <- function() {
  tag <- "gilad"
  tag <- "katzmouse"
  norm.norm <- c("tmm", "tmm", "deseq", "tmm")
  norm.test <- c("edger", "bayseq", "deseq", "deseq")
  test.test <- c("edger", "deseq", "bayseq")
  cols <- c("black", "red", "lightblue", "blue")
  labs <- c("iDEGES/edgeR", "iDEGES/TbT", "iDEGES/DESeq", "iDEGES/TDT")
  #lims <- list(c(0.005, 0.015), c(0, 0.01), c(0, 0.01))
  lims <- list(c(0.0, 0.015), c(0, 0.005), c(0, 0.005))

  for (t in 1:length(test.test)) {
    png(paste(tag, test.test[t], ".png", sep = "_"), 350, 280)
    par(mar = c(4, 4, 2, 2))
    plot(0, 0, type = "n", xlim = c(1, 10), ylim = lims[[t]], xlab = "iteration", ylab = "1 - correlation")
    for (n in 1:length(norm.norm)) {
      ranks <- read.table(paste("result_revise/", tag, "_RankVec_", norm.norm[n], "_",
                          norm.test[n], "_", test.test[t], "_.txt",sep = ""))
      y <- rep(0, length = ncol(ranks))
      for (r in 1:ncol(ranks)) {
        y[r] <- 1 - cor(ranks[, 1], ranks[, r])
      }
      lines(1:10, y[2:11], col = cols[n])
      points(1:10, y[2:11], col = cols[n])
    }
    #legend("topright", ncol = 2,
    #     legend = c("iDEGES/TbT", "iDEGES/edgeR", "iDEGES/TDT", "iDEGES/DESeq"),
    #     col = c("red", "black", "blue", "lightblue"), pch = 1, lty = 1, cex = 0.8)
    dev.off()
  }



}

clustalTasVector <- function() {
  tag <- "gilad"
  tag <- "katzmouse"
  count <- read.table(paste(data.path, "/", tag, "_count_table.txt", sep = ""), header = TRUE)
  rownames(count) <- count[, 1]
  count <- count[, - 1]
  count <- count[(rowSums(count) > 0), ]
  fmt <- read.table(paste("dentList_", tag, "2.txt", sep = ""), header = FALSE)
  fmt <- unique(fmt)
  df <- matrix(0, ncol = nrow(fmt), nrow = nrow(count))
  colnames(df) <- fmt[, 1]

  for (i in 1:nrow(fmt)) {
    f <- paste("result_revise/", tag,"_RankVec_", fmt[i, 2], "_", fmt[i, 3], "_", fmt[i, 5], "_.txt", sep = "")
    d <- read.table(f)
    e <- try(df[, i] <- d[, fmt[i, 4] + 1])
    if (class(e) == "try-error")
      browser()
  }

  mt <- matrix(0, ncol = nrow(fmt), nrow = nrow(fmt))
  colnames(mt) <- colnames(df)
  rownames(mt) <- colnames(df)
  for (i in 1:nrow(fmt)) {
    for (j in 1:nrow(fmt)) {
      if (i > j) {
        mt[i, j] <- 1 - cor(df[, i], df[, j])
      } else {
        mt[i, j] <- 0
      }
    }
  }
  d <- as.dist(mt)
  h <- hclust(d, method = "average")
  png(paste("clust", tag, "1.png", sep = "_"), width = 1040, height = 720)
  ma <- par(mar = c(1, 4, 1, 2))
  plot(h, ylab = "", xlab = "", main = "", sub = "", col = 1:10)
  dev.off()
}












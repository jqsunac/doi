# analyzeScript.R
#
# This script is for analyzing simulation results.
# Run this script after simulateScript.R can creates the part for Tables in paper.
# Usage:
#   R --slave --vanilla < analyzeScript.R
#
# It needs to set the directory which contains simulation results.
# The analyzed results are save to the same directory with prefix "anResult_".
result_path <- "./result"






analyzeAUC <- function(datadir, rep.type=NULL) {
  # analyzeAUC
  #
  # 1) Read simulation result(AUC) in "datadir" directory.
  # 2) Calculate the mean of AUC values of each conditions.
  # 3) Save the mean of AUC to file.
  # 4) Do t-paired test for all candidate situations.
  # 5) Adjust data and write to Table 1.
  #

  current.dir <- getwd()
  setwd(datadir)

  # Prepare the arguments for analyzing.
  files <- list.files()
  isAUC <- grep("auc", files)
  issimResult <- grep("simData_", files)
  isTarget <- grep(rep.type, files)
  files <- files[intersect(intersect(isAUC,issimResult),isTarget)]
  pdeg.label <- NULL
  pa.label <- NULL
  for (i in 1:length(files)) {
    rds <- strsplit(files[i], "_")
    pdeg.label <- union(pdeg.label, as.numeric(rds[[1]][5]))
    pa.label <- union(pa.label, as.numeric(rds[[1]][7]))
  }
  pdeg.label <- pdeg.label[!is.na(pdeg.label)]
  pa.label <- pa.label[!is.na(pa.label)]
  auc <- list()

  # Read all raw data and save it to auc list object.
  for (i in 1:length(files)) {
    d <- read.table(files[i], sep="\t", header=TRUE)
    rds <- strsplit(files[i], "_")
    pdeg.idx <- as.logical(rds[[1]][5] == pdeg.label)
    pa.idx <- as.logical(rds[[1]][7] == pa.label)
    methods <- colnames(d)
    for (j in 1:length(methods)) {
      if (is.null(auc[methods[j]][[1]]))
        auc[methods[j]][[1]] <- array(0, dim=(c(length(pdeg.label), length(pa.label), nrow(d))))
      auc[methods[j]][[1]][pdeg.idx, pa.idx, ] <- d[, j]
    }
  }
  # Compute the mean of AUC values and write to files.
  for (i in 1:length(auc)) {
    this.auc.mean <- matrix(0, nrow=length(pdeg.label), ncol=length(pa.label))
    rownames(this.auc.mean) <- pdeg.label
    colnames(this.auc.mean) <- pa.label
    this.auc.var <- matrix(0, nrow=length(pdeg.label), ncol=length(pa.label))
    rownames(this.auc.var) <- pdeg.label
    colnames(this.auc.var) <- pa.label
    for (j in 1:nrow(this.auc.mean)) {
      this.auc.mean[j, ] <- apply(auc[[i]][j, , ], 1, mean)
      this.auc.var[j, ] <- apply(auc[[i]][j, , ], 1, var)
    }
    write.table(this.auc.mean, file=paste("anResult_", rep.type, "_AUC_mean_", names(auc[i]), ".txt", sep=""),
      quote=FALSE, sep="\t", row.names=TRUE, col.name=TRUE)
    write.table(this.auc.var, file=paste("anResult_", rep.type, "_AUC_var_", names(auc[i]), ".txt", sep=""),
      quote=FALSE, sep="\t", row.names=TRUE, col.name=TRUE)
  }

  # Test the AUC values with wilcox.test two-side.
  for (i in 1:length(auc)) {
    for (j in i:length(auc)) {
      if (i != j) {
        test.result <- matrix("", nrow=length(pdeg.label), ncol=length(pa.label))
        rownames(test.result) <- pdeg.label
        colnames(test.result) <- pa.label
        for (pa.idx in 1:length(pa.label)) {
          for (pdeg.idx in 1:length(pdeg.label)) {
            r <- t.test(auc[[i]][pdeg.idx, pa.idx, ], auc[[j]][pdeg.idx, pa.idx, ], paired=TRUE)
            if (r$p.value < 0.01) {
              test.result[pdeg.idx, pa.idx] <- "+"
            } else {
              test.result[pdeg.idx, pa.idx] <- "-"
            }
          }
        }
        write.table(test.result, 
          file=paste("anResult_", rep.type, "_AUC_test_", names(auc[i]), "_vs_", names(auc[j]), ".txt", sep=""),
          quote=FALSE, sep="\t", row.names=TRUE, col.name=TRUE)
      }
    }
  }

  setwd(current.dir)
}

analyzeTime <- function(datadir, rep.type=NULL) {
  # analyzeTime
  #
  # 1) Read simulation result(execution time) in "datadir" directory.
  # 2) Calculate the rate of execution time.
  # 3) Get max and min values and write to files.
  #

  current.dir <- getwd()
  setwd(datadir)

  # Prepare the arguments for analyzing.
  files <- list.files()
  istime <- grep("time", files)
  issimResult <- grep("simData_", files)
  isTarget <- grep(rep.type, files)
  files <- files[intersect(intersect(istime,issimResult),isTarget)]
  pdeg.label <- NULL
  pa.label <- NULL
  for (i in 1:length(files)) {
    rds <- strsplit(files[i], "_")
    pdeg.label <- union(pdeg.label, as.numeric(rds[[1]][5]))
    pa.label <- union(pa.label, as.numeric(rds[[1]][7]))
  }
  extime <- list()
  
  # Read all raw data and save it to extime list object.
  times <- 0 # simulation times
  for (i in 1:length(files)) {
    d <- read.table(files[i], sep="\t", header=TRUE)
    if (i == 1)
      times <- nrow(d)
    rds <- strsplit(files[i], "_")
    pdeg.idx <- as.logical(rds[[1]][5] == pdeg.label)
    pa.idx <- as.logical(rds[[1]][7] == pa.label)
    methods <- colnames(d)
    for (j in 1:length(methods)) {
      if (is.null(extime[methods[j]][[1]]))
        extime[methods[j]][[1]] <- array(0, dim=(c(length(pdeg.label), length(pa.label), times)))
      extime[methods[j]][[1]][pdeg.idx, pa.idx, ] <- d[, j]
    }
  }

  # Compute the mean of AUC values and write to files.
  for (i in 1:length(extime)) {
    this.extime <- matrix(0, nrow=length(pdeg.label), ncol=length(pa.label))
    rownames(this.extime) <- pdeg.label
    colnames(this.extime) <- pa.label
    for (j in 1:nrow(this.extime)) {
      this.extime[j, ] <- apply(extime[[i]][j, , ], 1, mean)
    }
    write.table(this.extime, file=paste("anResult_", rep.type, "_time_mean_", names(extime[i]), ".txt", sep=""),
      quote=FALSE, sep="\t", row.names=TRUE, col.name=TRUE)
  }

  # Get the ratio between two methods.
  for (i in 1:length(extime)) {
    for (j in i:length(extime)) {
      if (i != j) {
        max.ratio <- matrix("", nrow=length(pdeg.label), ncol=length(pa.label))
        min.ratio <- matrix("", nrow=length(pdeg.label), ncol=length(pa.label))
        rownames(min.ratio) <- pdeg.label
        colnames(min.ratio) <- pa.label
        rownames(max.ratio) <- pdeg.label
        colnames(max.ratio) <- pa.label
        for (pa.idx in 1:length(pa.label)) {
          for (pdeg.idx in 1:length(pdeg.label)) {
            r <- extime[[i]][pdeg.idx, pa.idx, ] / extime[[j]][pdeg.idx, pa.idx, ]
            max.ratio[pdeg.idx, pa.idx] <- max(r)
            min.ratio[pdeg.idx, pa.idx] <- min(r)
          }
        }
        cat(paste("MAX (", names(extime[i]), "/", names(extime[j]), ") = ", max(max.ratio), "\n", sep=""))
        cat(paste("min (", names(extime[i]), "/", names(extime[j]), ") = ", max(min.ratio), "\n", sep=""))
      }
    }
  }

  setwd(current.dir)
}


analyzeEliminateACC <- function(datadir, rep.type=NULL) {
  # Analyze the accuracy of elimnation of step 2 of DES process.
  #
  #
  #
  current.dir <- getwd()
  setwd(datadir)
  
  # Prepare the arguments for analyzing.
  files <- list.files()
  issimResult <- grep("simData_", files)
  isTarget <- grep(rep.type, files)
  accfiles <- list()
  accfiles$TP <- files[intersect(intersect(issimResult, isTarget), grep("TP", files))]
  accfiles$FN <- files[intersect(intersect(issimResult, isTarget), grep("FN", files))]
  accfiles$FP <- files[intersect(intersect(issimResult, isTarget), grep("FP", files))]
  accfiles$TN <- files[intersect(intersect(issimResult, isTarget), grep("TN", files))]
  pdeg.label <- NULL
  pa.label <- NULL
  for (i in 1:length(accfiles$TP)) {
    rds <- strsplit(accfiles$TP[i], "_")
    pdeg.label <- union(pdeg.label, as.numeric(rds[[1]][5]))
    pa.label <- union(pa.label, as.numeric(rds[[1]][7]))
  }
  acc <- list()


  # Read data.
  for (i in 1:length(accfiles)) {
    acc[names(accfiles)[i]][[1]] <- list()

    for (j in 1:length(accfiles[[i]])) {
      d <- read.table(accfiles[[i]][j], header=T, sep="\t")
      rds <- strsplit(accfiles[[i]][j], "_")
      pdeg.idx <- as.logical(rds[[1]][5] == pdeg.label)
      pa.idx <- as.logical(rds[[1]][7] == pa.label)
      methods <- colnames(d)
      for (k in 1:length(methods)) {
        if (is.null(acc[names(accfiles)[i]][[1]][methods[k]][[1]]))
          acc[names(accfiles[i])][[1]][methods[k]][[1]] <- array(0, dim=c(length(pdeg.label), length(pa.label), nrow(d)))
        acc[names(accfiles[i])][[1]][methods[k]][[1]][pdeg.idx, pa.idx, ] <- d[, k]
      }
    }
  }

  # Write the average of TP, TN, FP, FN to files.
  for (i in 1:length(acc)) {
    methods <- names(acc[[i]])
    for (j in 1:length(methods)) {
      if (methods[j] == "TMM")
        next
      this.acc <- matrix(0, nrow=length(pdeg.label), ncol=length(pa.label))
      rownames(this.acc) <- pdeg.label
      colnames(this.acc) <- pa.label
      for (k in 1:nrow(this.acc)) {
        this.acc[k, ] <- apply(acc[[i]][methods[j]][[1]][k, , ], 1, mean)
      }
      # Write to file.
      # File name ex) "anResult_monorep_TN_mean_DES.txt"
      write.table(this.acc,
        file=paste("anResult_", rep.type, "_", names(acc[i]), "_mean_", methods[j], "_DESmiddle.txt", sep=""), 
        sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
    }
  }




  # Write the average ACC, TPR, FPR to files.
  for (i in 1:length(acc)) {
    methods <- names(acc[[i]])
    for (j in 1:length(methods)) {
      if (methods[j] == "TMM")
        next
      tmpl <- matrix(0, nrow=length(pdeg.label), ncol=length(pa.label))
      rownames(tmpl) <- pdeg.label
      colnames(tmpl) <- pa.label
      this.acc <- tmpl
      this.tpr <- tmpl
      this.fpr <- tmpl
      for (k in 1:nrow(this.acc)) {
        this.acc[k, ] <- apply(
            (acc$TP[methods[j]][[1]][k, , ] + acc$TN[methods[j]][[1]][k, , ]),
            1, mean)
        this.tpr[k, ] <- apply(
            (acc$TP[methods[j]][[1]][k, , ] / (acc$TP[methods[j]][[1]][k, , ] + acc$FN[methods[j]][[1]][k, , ])),
            1, mean)
        this.fpr[k, ] <- apply(
            (acc$FP[methods[j]][[1]][k, , ] / (acc$FP[methods[j]][[1]][k, , ] + acc$TN[methods[j]][[1]][k, , ])),
            1, mean)
      }
      write.table(this.acc,
        file=paste("anResult_", rep.type, "_ACC_mean_", methods[j], "_DESmiddle.txt", sep=""), 
        sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
      write.table(this.tpr,
        file=paste("anResult_", rep.type, "_TPR_mean_", methods[j], "_DESmiddle.txt", sep=""), 
        sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
      write.table(this.fpr,
        file=paste("anResult_", rep.type, "_FPR_mean_", methods[j], "_DESmiddle.txt", sep=""), 
        sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
    }
  }
  setwd(current.dir)

}



cat("--------MULTI REPLICATES----------------------\n")
analyzeAUC(result_path, rep.type="multirep")
analyzeTime(result_path, rep.type="multirep")
analyzeEliminateACC(result_path, rep.type="multirep")
cat("--------MONO REPLICATE----------------------\n")
analyzeAUC(result_path, rep.type="monorep")
analyzeTime(result_path, rep.type="monorep")
analyzeEliminateACC(result_path, rep.type="monorep")




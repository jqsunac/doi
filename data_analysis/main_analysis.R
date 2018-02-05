## R code for homeolog expression analysis

set.seed(2016)
suppressMessages(require(MASS))
#suppressMessages(require(RUnit))
#suppressMessages(require(knitr))
suppressMessages(require(XLConnect))
suppressMessages(require(RColorBrewer))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))
suppressMessages(require(ggsci))
suppressMessages(require(ggExtra))
suppressMessages(require(gplots))
#suppressMessages(require(VennDiagram))
#suppressMessages(require(grid))
#suppressMessages(require(scales))
#suppressMessages(require(outliers))
#suppressMessages(require(fBasics))
#suppressMessages(require(edgeR))
#suppressMessages(require(TCC))
suppressMessages(require(GO.db))
suppressMessages(require(org.At.tair.db))
suppressMessages(require(clusterProfiler))
#suppressMessages(require(plyr))
#suppressMessages(require(genefilter))
#suppressMessages(require(googleVis))
#suppressMessages(require(venneuler))


options(stringsAsFactors = FALSE)


## GLOBAL PARAMTERS
FPKM_CUTOFF <- 1.0


## INIT DIRECTORIES
WORK_DIR <- dirname('./')
PATH <- list(
    res = paste0(WORK_DIR, "/result_files"),
    lib = paste0(WORK_DIR, "/lib"),
    src = paste0(WORK_DIR, "/lib/src"),
    dat = paste0(WORK_DIR, "/data"),
    tmp = paste0(WORK_DIR, "/tmp")
)
for (.p in PATH) if(!file.exists(.p)) dir.create(.p)
path.lib      <- paste0(PATH$res, "/library_features")
path.de       <- paste0(PATH$res, "/differential_expression")
path.geneprof <- paste0(PATH$res, "/expression_profile")
path.expbias  <- paste0(PATH$res, "/expression_bias")
path.expbias.path      <- paste0(path.expbias, "/path")
for (.path.i in c(path.lib, path.de, path.geneprof, path.expbias, path.expbias.path)) {
    if (!file.exists(.path.i)) dir.create(.path.i)
}







## INIT COLOR SETTINGS
.brewer.set1    <- brewer.pal(9, "Set1")
.brewer.pastel1 <- brewer.pal(9, "Pastel1")
.brewer.dark2   <- brewer.pal(8, "Dark2")
.brewer.greens  <- brewer.pal(9, "Greens")
COLS <- list(
    ama  = .brewer.set1[2],     # C. amara
    riv  = .brewer.set1[1],     # C. rivularis
    hir  = .brewer.set1[1],     # C. hirsuta
    ins  = .brewer.set1[3],     # C. insueta
    insA = .brewer.set1[4],     # C. insueta (A)
    insR = .brewer.set1[5],     # C. insueta (R)
    DE   = .brewer.set1[5],     # Differentially expressed genes
    NDE  = .brewer.set1[9]      # Non-differentially expressed genes
)

COLSFUNC <- function(n = 100) {
    #.colfunc <- colorRampPalette(c("#FFFFFF", "#D6604D", "#67000D"))
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu")))
    .colfunc(n)
}
COLSFUNCSP <- function(n = 100) {
    #.colfunc <- colorRampPalette(c("#FFFFFF", "#D6604D", "#67000D"))
    .colfunc <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
    .colfunc(n)
}
COLSFUNC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(c("#053061", "#4393C3", "#FFFFFF", "#D6604D", "#67001F"))
    .colfunc(n)
}
COLSFUNCC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu"))[5:9])
    .colfunc(n)
}



SPLABELS  <- c("AA", "IA", "IR", "RR")
SPLABELS2 <- c("AA", "IA", "II", "IR", "RR")
TIMELABELS  <- c("0 hr", "2 hr", "4 hr", "8 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr")
TIMELABELS2 <-  c("0 hr vs. 2 hr",   "0 hr vs. 4 hr",  "0 hr vs. 8 hr", "0 hr vs. 12 hr",
                  "0 hr vs. 24 hr", "0 hr vs. 48 hr", "0 hr vs. 72 hr", "0 hr vs. 96 hr")
INCH2CM <- 1 / 2.54





#' Split string of TAIR ID and convert them to CARHR ID.
#'
#'
TAIRSTR2CARHRVEC <- function(x, sep = '/') {
    paste0(unlist(TAIR2CARHR[unlist(strsplit(x, '/'))]), collapse = sep)
}



#' Plot multiple ggplot objects in a figure
#'
#'
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#' Save list object into Excel file.
#'
#'
SaveExcel <- function(data.list, file.name = NULL) {
    dat <- list()
    nam <- NULL
    for (i in names(data.list)) {
        if (class(data.list[[i]]) == "list") {
            for (j in names(data.list[[i]])) {
                nam <- c(nam, paste0(i, ".", j))
                dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]][[j]]), data.list[[i]][[j]])))
            }
        } else {
            nam <- c(nam, i)
            dat <- c(dat, list(data.frame(Rowname = rownames(data.list[[i]]), data.list[[i]])))
        }
    }
    names(dat) <- nam
    if (is.null(file.name)) stop("Need file name.")
    if (is.null(names(dat))) stop("Need list names.")
    if (file.exists(file.name)) file.remove(file.name)
    gc()
    workbook <- loadWorkbook(file.name, create = TRUE)
    createSheet(workbook, names(dat))
    writeWorksheet(workbook, dat, names(dat), header = TRUE)
    saveWorkbook(workbook)
    gc()
}





GetDesign <- function(tissue = "leaf", species = NULL) {
    design <- read.table(paste0(PATH$lib, "/src/lib_design.txt"), header = TRUE)
    keeped.libs <- (design$tissue == tissue)
    design      <- design[keeped.libs, ]
    libs.order  <- order(design$tissue, design$species, design$time)
    design      <- design[libs.order, ]
    design.sp   <- NULL
    if (!is.null(species)) {
        for (sp in species) {
            design.sp <- rbind(design.sp, design[grep(sp, design$species), ])
        }
        design <- design.sp
        design$species <- factor(design$species, levels = species)
    }
    design
}




GetCounts <- function(ca, design) {
    counts.aorigin <- counts.rorigin <- counts.common <- 
        vector("list", length = length(design$library))
    names(counts.aorigin) <- names(counts.rorigin) <- names(counts.common) <-
        design$library
    
    gene.names <- NULL
    
    for (libname in design$library) {
        .a <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_orig.sorted.txt"),
                         sep = "\t", fill = TRUE)
        .r <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_other.sorted.txt"),
                         sep = "\t", fill = TRUE)
        .c <- read.table(paste0(PATH$dat, "/counts/", libname,
                                "__on__camara.filtered.rc_common.sorted.txt"),
                         sep = "\t", fill = TRUE)
        gene.names <- .a[, 1]
        counts.aorigin[[libname]] <- .a[, 2]
        counts.rorigin[[libname]] <- .r[, 2]
        counts.common[[libname]] <- .c[, 2]
    }
    
    counts.aorigin <- as.data.frame(counts.aorigin)[-c(grep("^__", gene.names)), ]
    counts.rorigin <- as.data.frame(counts.rorigin)[-c(grep("^__", gene.names)), ]
    counts.common  <- as.data.frame(counts.common)[-c(grep("^__", gene.names)), ]
    
    print(data.frame(A      = colSums(counts.aorigin),
                     common = colSums(counts.common),
                     R      = colSums(counts.rorigin)))
    
    gene.names  <- gene.names[-c(grep("^__", gene.names))]
    
    counts.AA <- counts.aorigin[, 1:9] + counts.rorigin[, 1:9] + counts.common[, 1:9]
    counts.RR <- counts.aorigin[, 19:27] + counts.rorigin[, 19:27] + counts.common[, 19:27]
    
    aorigin.ratio <- (counts.aorigin / (counts.aorigin + counts.rorigin))[, 10:18]
    aorigin.ratio[(counts.aorigin[, 10:18] + counts.rorigin[, 10:18]) == 0 & counts.common[, 10:18] != 0] <- 1/3
    counts.IA <- counts.aorigin[, 10:18] + counts.common[, 10:18] * aorigin.ratio
    counts.IR <- counts.rorigin[, 10:18] + counts.common[, 10:18] * (1 - aorigin.ratio)
    counts.IA[is.na(counts.IA)] <- 0
    counts.IR[is.na(counts.IR)] <- 0
    
    ## add gene names and time course labels
    colnames(counts.AA) <- colnames(counts.RR) <- colnames(counts.IA) <- colnames(counts.IR) <-
            colnames(aorigin.ratio) <- TIMELABELS
    rownames(counts.AA) <- rownames(counts.RR) <- rownames(counts.IA) <- rownames(counts.IR) <-
            rownames(aorigin.ratio) <- gene.names
    
    ## remove all zero rows
    is.allzero <- as.logical(rowSums(counts.aorigin + counts.rorigin + counts.common) == 0)
    gene.names <- gene.names[!is.allzero]
    counts.AA <- counts.AA[!is.allzero, ]
    counts.RR <- counts.RR[!is.allzero, ]
    counts.IA <- counts.IA[!is.allzero, ]
    counts.IR <- counts.IR[!is.allzero, ]
    aorigin.ratio <- aorigin.ratio[!is.allzero, ]
    
    ca$counts <- list(AA = as.matrix(counts.AA), IA = as.matrix(counts.IA),
                      IR = as.matrix(counts.IR), RR = as.matrix(counts.RR))
    ca$aorigin.ratio <- as.matrix(aorigin.ratio)
    ca$gene.names <- gene.names
    ca
}

CalcFPKM <- function(ca) {
    gene.length <- read.table(paste0(PATH$lib, "/src/gene_length.tsv"), header = FALSE)
    row.names(gene.length) <- gene.length[, 1]
    gene.length <- gene.length[ca$gene.names, -1]
    
    fpkm <- ca$counts
    
    ## AA
    cpm.A <- sweep(ca$counts$AA, 2, 1e6 / colSums(ca$counts$AA), "*")
    fpkm$AA <- sweep(cpm.A, 1, 1e3 / gene.length[, 1], "*")
    
    ## RR
    cpm.R <- sweep(ca$counts$RR, 2, 1e6 / colSums(ca$counts$RR), "*")
    fpkm$RR <- sweep(cpm.R, 1, 1e3 / gene.length[, 1], "*")
    
    ## IA, IR
    print("IA counts ---> IA FPKM")
    print("IR counts ---> IR FPKM")
    print("IA counts + IR counts---> II FPKM")
    cpm.IA <- sweep(ca$counts$IA, 2, 1e6 / colSums(ca$counts$IA), "*")
    fpkm$IA <- sweep(cpm.IA, 1, 1e3 / gene.length[, 1], "*")
    cpm.IR <- sweep(ca$counts$IR, 2, 1e6 / colSums(ca$counts$IR), "*")
    fpkm$IR <- sweep(cpm.IR, 1, 1e3 / gene.length[, 1], "*")
    .counts.I <- ca$counts$IA + ca$counts$IR
    cpm.II <- sweep(.counts.I, 2, 1e6 / colSums(.counts.I), "*")
    fpkm$II <- sweep(cpm.II, 1, 1e3 / gene.length[, 1], "*")
   
    ca$gene.length <- gene.length
    ca$fpkm <- fpkm
    ca
}





DoGO2 <- function(sig.genes = NULL, all.genes = NULL, ontology = "BP",
                  p.cutoff = 1, q.cutoff = 1) {
    all.genes <- as.character(all.genes)
    sig.genes <- as.character(sig.genes)
    all.genes.at <- as.character(unlist(CARHR2TAIR[all.genes]))
    sig.genes.at <- as.character(unlist(CARHR2TAIR[sig.genes]))
    all.genes.at <- all.genes.at[all.genes.at != ""]
    sig.genes.at <- sig.genes.at[sig.genes.at != ""]
    ego <- enrichGO(gene = sig.genes.at, universe = all.genes.at, OrgDb = org.At.tair.db,
                    ont = ontology, pAdjustMethod = "BH", keyType = "TAIR",
                    pvalueCutoff = p.cutoff, qvalueCutoff = q.cutoff, readable = FALSE)
    egos <- dropGO(ego, level = 1)
    egos <- dropGO(egos, level = 2)
    egos <- dropGO(egos, level = 3)
    egostable <- as.data.frame(egos)
    numgenes <- sapply(strsplit(egostable$BgRatio, "/"), function(x) as.integer(x[1]))
    egostable <- egostable[(10 < numgenes) & (numgenes < 500), ]
    egostable <- egostable[, colnames(egostable) != "p.adjust"]
    egostable$qvalue <- p.adjust(egostable$pvalue, method = "BH")
    
    .carhr <- rep("", length = nrow(egostable))
    for (.i in 1:nrow(egostable)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(egostable$geneID[.i])
    }
    egostable <- data.frame(egostable, carID = .carhr)
    
    list(GO = ego, GOTABLE = as.data.frame(ego), GOSMPLTABLE = egostable)
}









LargeFPKM_GO <- function(ca) {
    A <- ca$fpkm$AA
    R <- ca$fpkm$RR
    I <- ca$fpkm$IA + ca$fpkm$IR
    A <- A[rowSums(A > FPKM_CUTOFF) > 0, ]
    R <- R[rowSums(R > FPKM_CUTOFF) > 0, ]
    I <- I[rowSums(I > FPKM_CUTOFF) > 0, ]
    A <- log10(A+1)
    R <- log10(R+1)
    I <- log10(I+1)
    sd.A <- apply(A, 1, sd, na.rm = T)
    sd.R <- apply(R, 1, sd, na.rm = T)
    sd.I <- apply(I, 1, sd, na.rm = T)
    mu.A <- apply(A, 1, mean, na.rm = T)
    mu.R <- apply(R, 1, mean, na.rm = T)
    mu.I <- apply(I, 1, mean, na.rm = T)
    cv.A <- sd.A / mu.A
    cv.R <- sd.R / mu.R
    cv.I <- sd.I / mu.I
    
    # camara
    kamara <- cv.A > 0.20 & mu.A > 1
    goobj <- DoGO(sig.genes = names(kamara[kamara]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableama<- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableama))
    for (.i in 1:nrow(gotableama)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableama$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableama, .carhr)),
              file.name = paste0(path.expbias, "/GO_camara_CV0.20_and_mu1.0.xlsx"))
    gc()
    # crivularis
    krivularis <- cv.R > 0.20 & mu.R > 1
    goobj <- DoGO(sig.genes = names(krivularis[krivularis]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableriv <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableriv))
    for (.i in 1:nrow(gotableriv)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableriv$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableriv, .carhr)),
              file.name = paste0(path.expbias, "/GO_rivularis_CV0.20_and_mu1.0.xlsx"))
    gc()
    # cinsueta
    kinsueta <- cv.I > 0.20 & mu.I > 1
    goobj <- DoGO(sig.genes = names(kinsueta[kinsueta]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableins <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableins))
    for (.i in 1:nrow(gotableins)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableins$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableins,.carhr)),
              file.name = paste0(path.expbias, "/GO_cinuseta_CV0.20_and_mu1.0.xlsx"))
    gc()
    
    # cinsueta (var AOriginRatio)
    theta.ar <- ca$aorigin.ratio[names(mu.I), ]
    var.theta.ar <- apply(theta.ar, 1, var, na.rm = T)
    keep2c <- var.theta.ar > 0.01  & mu.I > 1
    goobj <- DoGO(sig.genes = names(keep2c[keep2c]), all.genes = ca$EXP$genes$ALL, p.cutoff = 1, q.cutoff = 1)
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    .target.terms <- intersect(goobj$TABLE$ID, target.terms)
    gotablesimple <- goobj$TABLE[.target.terms,]
    gotablesimple$qvalue <- p.adjust(gotablesimple$pvalue, method = "BH")
    gotableins2 <- gotablesimple[, -6]
    .carhr <- rep("", length=nrow(gotableins2))
    for (.i in 1:nrow(gotableins2)) {
        .carhr[.i] <- TAIRSTR2CARHRVEC(gotableins2$geneID[.i])
    }
    SaveExcel(list(full=goobj$TABLE, simlify=cbind(gotableins2, .carhr)),
              file.name = paste0(path.expbias, "/GO_cinsueta_VarARatio0.01_and_mu1.0.v2.xlsx"))

    
    goterms <- unique(c(gotableama$ID[gotableama$qvalue < 0.1],
                        gotableriv$ID[gotableriv$qvalue < 0.1],
                        gotableins$ID[gotableins$qvalue < 0.1],
                        gotableins2$ID[gotableins2$qvalue < 0.1]))
    
    go2desc <- rep("", length=length(goterms))
    gomat <- matrix(NA, ncol = 4, nrow = length(goterms))
    rownames(gomat) <- names(go2desc) <- goterms
    colnames(gomat) <- c("AA_CvMean", "II_CvMean", "II_VRatioMean", "RR_CvMean")
    go2desc[gotableama$ID[gotableama$qvalue < 0.1]] <- gotableama$Description[gotableama$qvalue < 0.1]
    go2desc[gotableriv$ID[gotableriv$qvalue < 0.1]] <- gotableriv$Description[gotableriv$qvalue < 0.1]
    go2desc[gotableins$ID[gotableins$qvalue < 0.1]] <- gotableins$Description[gotableins$qvalue < 0.1]
    go2desc[gotableins2$ID[gotableins2$qvalue < 0.1]] <- gotableins2$Description[gotableins2$qvalue < 0.1]
    
    gomat[gotableama$ID[gotableama$qvalue < 0.1], 1] <- -log10(gotableama$qvalue)[gotableama$qvalue < 0.1]
    gomat[gotableins$ID[gotableins$qvalue < 0.1], 2] <- -log10(gotableins$qvalue)[gotableins$qvalue < 0.1]
    gomat[gotableins2$ID[gotableins2$qvalue < 0.1], 3] <- -log10(gotableins2$qvalue)[gotableins2$qvalue < 0.1]
    gomat[gotableriv$ID[gotableriv$qvalue < 0.1], 4] <- -log10(gotableriv$qvalue)[gotableriv$qvalue < 0.1]
    
    
    gomat2 <- gomat
    gomat2[is.na(gomat2)] <- 1
    h <- heatmap.2(gomat2, trace = "none", scale = "none", dendrogram = "row", Colv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNCC2(), key.xlab = "log2(FPKM + 1)")
    df <- data.frame(gomat, desc = go2desc)
    df <- df[h$rowInd, ]
    SaveExcel(list(Sheet1=df), file.name = paste0(path.lib, "/CVLogFPKM0.2_meanLogFPKM1.0.xlsx"))
}






PlotRqsMDS <- function(ca) {
    rsqmat <- matrix(NA, ncol = 4 * 9, nrow = 4 * 9)
    f <- cbind(ca$fpkm$AA, ca$fpkm$IA, ca$fpkm$IR, ca$fpkm$RR)
    f <- log10(f[ca$EXP$gene$ALL, ] + 1)
    
    sd.log10.f <- apply(f, 1, sd, na.rm = T)
    mu.log10.f <- apply(f, 1, mean, na.rm = T)
    cv.log10.f <- sd.log10.f / mu.log10.f
    
    #f <- f[(cv.log10.f > 1.5), ]
    
    for (i in 2:(4*9)) {
        for (j in 1:(i-1)) {
            rsqmat[i, j] <- summary(lm(f[, i] ~ f[, j]))$r.squared
        }
    }
    
    colnames(rsqmat) <- rownames(rsqmat) <- 
        paste0(rep(SPLABELS, each = 9), "__", rep(TIMELABELS, times = 4))
    
    mds <- data.frame(cmdscale(as.dist(1 - rsqmat)))
    mdsdf <- data.frame(
        x = mds$X1, y = mds$X2,
        species = sapply(strsplit(rownames(mds), "__"), function(x) x[1]),
        time = gsub(" hr", "", sapply(strsplit(rownames(mds), "__"), function(x) x[2]))
    )
    g <- ggplot(mdsdf, aes(x = x, y = y, color = species, label = time))
    g <- g + geom_text(size=6) + coord_fixed()
    g <- g + xlab("X1") + ylab("X2")
    g <- g + scale_color_manual(values = c(COLS$ama, COLS$insA, COLS$insR, COLS$riv))
    g <- g + theme_bw()
    g <- g + xlim(-0.5, 0.5) + ylim(-0.5, 0.5)
    g
    
    png(paste0(path.lib, "/logFPKM.MDS.png"), 340, 340)
    plot(g)
    dev.off()
}


PlotPCA <- function(ca) {

    f <- cbind(ca$fpkm$AA, ca$fpkm$IA, ca$fpkm$IR, ca$fpkm$RR)
    f <- log10(f[ca$EXP$gene$ALL, ] + 1)
    
    fpkm.class <- rep(c('A-FPKM', 'IA-FPKM', 'IR-FPKM', 'R-FPKM'), each = 9)
    time <- gsub(' hr', '', colnames(f))
    colnames(f) <- paste0(fpkm.class, ' ', time)

    pcaobj <- prcomp(t(f), scale = FALSE)
    # plot proportion of variances of PCA
    pca.prop.vars <- data.frame(Var = summary(pcaobj)$importance[2, ],
                                PC = factor(colnames(summary(pcaobj)$importance),
                                            levels = colnames(summary(pcaobj)$importance)))
    pca.prop.vars.ggplot <- ggplot(pca.prop.vars, aes(x = PC, y = Var))
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + geom_bar(stat = 'identity')
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + theme_bw()
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + xlab('principal component')
    pca.prop.vars.ggplot <- pca.prop.vars.ggplot + ylab('proportion of variance')
    png(paste0(path.lib, '/pca.prop.vars.png'), 800, 400)
    plot(pca.prop.vars.ggplot)
    dev.off()
    
    pc.coord <- t(t(pcaobj$rotation) %*% f)
    pca.pc.coords <- data.frame(pc.coord,
                                fpkm = fpkm.class,
                                time = time)
    pca.pc.coords.ggplot <- ggplot(pca.pc.coords, aes(x = PC1, y = PC2, color = fpkm, label = time))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + geom_text() + scale_color_nejm() + coord_fixed()
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + theme_bw() + guides(color = guide_legend(title = ''))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlab(paste0('PC1 (', summary(pcaobj)$importance[2, 1] * 100, '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + ylab(paste0('PC2 (', summary(pcaobj)$importance[2, 2] * 100, '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlim(-20, 40) + ylim(-70, -10)
    png(paste0(path.lib, '/pca.pc1pc2.png'), 400, 400)
    plot(pca.pc.coords.ggplot)
    dev.off()

    pca.pc.coords.ggplot <- ggplot(pca.pc.coords, aes(x = PC2, y = PC3, color = fpkm, label = time))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + geom_text() + scale_color_nejm() + coord_fixed()
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + theme_bw() + guides(color = guide_legend(title = ''))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlab(paste0('PC2 (', summary(pcaobj)$importance[2, 2] * 100, '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + ylab(paste0('PC3 (', summary(pcaobj)$importance[2, 3] * 100, '%)'))
    pca.pc.coords.ggplot <- pca.pc.coords.ggplot + xlim(-70, -10) + ylim(-40, 20)
    png(paste0(path.lib, '/pca.pc2pc3.png'), 400, 400)
    plot(pca.pc.coords.ggplot)
    dev.off()

}





PlotAoriginRatioDist <- function(ca) {
    r <- ca$aorigin.ratio
    target.homeologs <- (ca$fpkm$IA + ca$fpkm$IR > FPKM_CUTOFF)
    dfr <- NULL
    for (i in 1:9) {
        ri <- r[target.homeologs[, i], i]
        dfi <- data.frame(value = ri, gene = names(ri), time = paste0(TIMELABELS[i], " (", length(ri), " homeologs)"))
        dfr <- rbind(dfr, dfi)
    }
    
    dfr$time <- factor(dfr$time, levels = paste0(TIMELABELS, " (", apply(target.homeologs, 2, sum), " homeologs)"))
    
    g <- ggplot(dfr, aes(x = value))
    g <- g + geom_histogram(binwidth = 0.025)
    g <- g + theme_bw()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "white"))
    g <- g + scale_x_continuous(breaks = c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1), minor_breaks = NULL)
    g <- g + facet_grid(~ time)#, ncol = 5)
    g <- g + xlab("A-origin ratio") + ylab("frequency")
    ca$BIAS$aorigin.ratio.dist <- g

    png(paste0(path.expbias, "/AOriginRatio.hist.lib.png"), 1000, 240)
    options(digits=2)
    plot(g)
    dev.off()
    options(digits=6)
     
     
    dfr <- NULL
    for (i in 1:9) {
        ri <- r[target.homeologs[, i], i]
        dfi <- data.frame(value = ri, gene = names(ri), time = TIMELABELS[i])
        dfr <- rbind(dfr, dfi)
    }
    
    gene2chr_txt <- read.table(paste0(PATH$lib, "/src/gene_chr.txt"), header = F, sep = "\t")
    gene2chr <- gene2chr_txt[, 2]
    names(gene2chr) <- gene2chr_txt[, 1]
    dfr <- data.frame(dfr, chr = gene2chr[dfr$gene])
    dfr <- dfr[grep("Chr", dfr$chr), ]
    g <- ggplot(dfr, aes(x = value))
    g <- g + geom_histogram(binwidth = 0.025)
    g <- g + theme_bw()
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "white"))
    g <- g + scale_x_continuous(breaks = c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1), minor_breaks = NULL)
    g <- g + facet_grid(chr ~ time, scale="free_y")
    g <- g + xlab("A-origin ratio") + ylab("frequency")
    png(paste0(path.expbias, "/AOriginRatio.hist.chr.png"), 1500, 1200)
    plot(g)
    dev.off()
    
    ca
}





PlotProfile <- function(g, ca, cls = TRUE, m = c(10, 10), col = NULL) {
    g <- intersect(g, rownames(ca$fpkm$AA))
    g.fpkm.A <- ca$fpkm$AA[g, ]
    g.fpkm.a <- ca$fpkm$IA[g, ]
    g.fpkm.r <- ca$fpkm$IR[g, ]
    g.fpkm.R <- ca$fpkm$RR[g, ]
    g.fpkm <- cbind(cbind(g.fpkm.A, g.fpkm.a), cbind(g.fpkm.r, g.fpkm.R))
    df <- BindTairName(g)[, -4]
    df[df[, 3] == "NA", 3] <- df[df[, 3] == "NA", 2]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 2]
    df[, 3] <- gsub("\\|\\|", ",", df[, 3])
    df[df[, 2] == "", 2] <- df[df[, 2] == "", 1]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 1]
    if (!is.null(col)) df[, 3] <- paste0("[", col, "] ", df[, 3])
    colnames(g.fpkm) <- paste0(rep(SPLABELS, each = 9), " (", TIMELABELS, ")")
    rownames(g.fpkm) <- df[, 3]
    g.fpkm <- log10(g.fpkm + 1)
    if (cls) {
        heatmap.2(g.fpkm, trace = "none", scale = "row", dendrogram = "row", Colv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "z-score", margins = m,
              colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "z-score", margins = m,
              colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    }
}


PlotProfile.FC <- function(g, ca, cls = TRUE, m = c(10, 10), col = NULL) {
    g <- intersect(g, rownames(ca$fpkm$AA))
    g.fpkm.A <- log2(ca$fpkm$AA[g, ] + 1)
    g.fpkm.a <- log2(ca$fpkm$IA[g, ] + 1)
    g.fpkm.r <- log2(ca$fpkm$IR[g, ] + 1)
    g.fpkm.R <- log2(ca$fpkm$RR[g, ] + 1)
    g.fpkm.A <- (g.fpkm.A - g.fpkm.A[, 1])#[, -1]
    g.fpkm.a <- (g.fpkm.a - g.fpkm.a[, 1])#[, -1]
    g.fpkm.r <- (g.fpkm.r - g.fpkm.r[, 1])#[, -1]
    g.fpkm.R <- (g.fpkm.R - g.fpkm.R[, 1])#[, -1]

    g.fpkm <- cbind(cbind(g.fpkm.A, g.fpkm.a), cbind(g.fpkm.r, g.fpkm.R))
    df <- BindTairName(g)[, -4]
    df[df[, 3] == "NA", 3] <- df[df[, 3] == "NA", 2]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 2]
    df[, 3] <- gsub("\\|\\|", ",", df[, 3])
    df[df[, 2] == "", 2] <- df[df[, 2] == "", 1]
    df[df[, 3] == "", 3] <- df[df[, 3] == "", 1]
    if (!is.null(col)) df[, 3] <- paste0("[", col, "] ", df[, 3])

    splabels <- c('C. amara', 'C. insueta (A-sub)', 'C. insueta (R-sub)', 'C. rivularis')
    colnames(g.fpkm) <- paste0(rep(splabels, each = 9), " (", TIMELABELS, ")")
    
    rownames(g.fpkm) <- df[, 3]
    g.fpkm <- g.fpkm[, -c(1, 9+1, 18+1, 27+1)]
    if (cls) {
        h <- heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "row", Colv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        h <- heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    }
    invisible(h)
}




CalcExpressedGenes <- function(ca) {
    exp.AA <- rownames(ca$fpkm$AA)[apply(ca$fpkm$AA > FPKM_CUTOFF, 1, any)]
    exp.RR <- rownames(ca$fpkm$RR)[apply(ca$fpkm$RR > FPKM_CUTOFF, 1, any)]
    exp.IA <- rownames(ca$fpkm$IA)[apply(ca$fpkm$IA > FPKM_CUTOFF, 1, any)]
    exp.IR <- rownames(ca$fpkm$IR)[apply(ca$fpkm$IR > FPKM_CUTOFF, 1, any)]
    exp.all <- unique(c(exp.AA, exp.RR, exp.IA, exp.IR))
    print(paste0("Overlap(AA, RR): ",  round(length(intersect(exp.AA, exp.RR))/length(union(exp.AA, exp.RR)), 3)))
    print(paste0("Overlap(IA, IR): ",  round(length(intersect(exp.IA, exp.IR))/length(union(exp.IA, exp.IR)), 3)))
    print(paste0("Overlap(AA, IA): ",  round(length(intersect(exp.AA, exp.IA))/length(union(exp.AA, exp.IA)), 3)))
    print(paste0("Overlap(RR, IR): ",  round(length(intersect(exp.RR, exp.IR))/length(union(exp.RR, exp.IR)), 3)))
    
    png(paste0(path.lib, "/Exp-genes-venn.png"), 500, 500)
    v <- venn(list(AA = exp.AA, IA = exp.IA, IR = exp.IR, RR = exp.RR))
    dev.off()
    
    genes <- gofull <- gosimple <- vector("list", length = length(attr(v, "intersections")))
    names(genes) <- names(gofull) <- names(gosimple) <- c(names(attr(v, "intersections")))
    for (vn in names(attr(v, "intersections"))) {
        .g <- attr(v, "intersections")[[vn]]
        .go <- DoGO2(sig.genes = .g, all.genes = exp.all)
        gofull[[vn]] <- .go$GOTABLE
        gosimple[[vn]] <- .go$GOSMPLTABLE
        genes[[vn]] <- BindTairName(.g)
    }
    
    names(gofull)    <- gsub(":", "-", names(gofull))
    names(gosimple) <- gsub(":", "-", names(gosimple))
    names(genes)     <- gsub(":", "-", names(genes))
    SaveExcel(gofull, paste0(path.lib, "/Exp-genes-GO-full.xlsx"))
    gc()
    SaveExcel(gosimple, paste0(path.lib, "/Exp-genes-GO-simple.xlsx"))
    gc()
    SaveExcel(genes, paste0(path.lib, "/Exp-genes.xlsx"))
    gc() 
    
    # expressed genes
    ca$EXP$genes <- list(AA = exp.AA, RR = exp.RR,
                         IA = exp.IA, IR = exp.IR,
                         ALL = exp.all)
    
    neo.exp.IA <- setdiff(exp.IA, exp.AA)
    neo.exp.IR <- setdiff(exp.IR, exp.RR)
    sup.exp.IA <- setdiff(exp.AA, exp.IA)
    sup.exp.IR <- setdiff(exp.RR, exp.IR)
    
    print(paste0("Newly expressed in A: ", length(neo.exp.IA)))
    print(paste0("Suppressed in A:      ", length(sup.exp.IA)))
    print(paste0("Newly expressed in R: ", length(neo.exp.IR)))
    print(paste0("Suppressed in R:      ", length(sup.exp.IR)))
    
    print(length(intersect(neo.exp.IA, neo.exp.IR))/length(union(neo.exp.IA, neo.exp.IR)))
    print(length(intersect(sup.exp.IA, sup.exp.IR))/length(union(sup.exp.IA, sup.exp.IR)))
    
    neo.exp.IA.go <- DoGO2(sig.genes = neo.exp.IA, all.genes = exp.all)
    neo.exp.IR.go <- DoGO2(sig.genes = neo.exp.IR, all.genes = exp.all)
    sup.exp.IA.go <- DoGO2(sig.genes = sup.exp.IA, all.genes = exp.all)
    sup.exp.IR.go <- DoGO2(sig.genes = sup.exp.IR, all.genes = exp.all)
    
    neo.IA.df <- list(gene = BindTairName(neo.exp.IA),
                      GOFULL = neo.exp.IA.go$GOTABLE,
                      GOSMPL = neo.exp.IA.go$GOSMPLTABLE)
    neo.IR.df <- list(gene = BindTairName(neo.exp.IR),
                      GOFULL = neo.exp.IR.go$GOTABLE,
                      GOSMPL = neo.exp.IR.go$GOSMPLTABLE)
    sup.IA.df <- list(gene = BindTairName(sup.exp.IA),
                      GOFULL= sup.exp.IA.go$GOTABLE,
                      GOSMPL = sup.exp.IA.go$GOSMPLTABLE)
    sup.IR.df <- list(gene = BindTairName(sup.exp.IR),
                      GOFULL= sup.exp.IR.go$GOTABLE,
                      GOSMPL = sup.exp.IR.go$GOSMPLTABLE)
    
    SaveExcel(neo.IA.df, paste0(path.lib, "/Newly-expressed-genes-A-genome.xlsx"))
    gc()
    SaveExcel(sup.IA.df, paste0(path.lib, "/Suppresed-genes-A-genome.xlsx"))
    gc()
    SaveExcel(neo.IR.df, paste0(path.lib, "/Newly-expressed-genes-R-genome.xlsx"))
    gc()
    SaveExcel(sup.IR.df, paste0(path.lib, "/Suppresed-genes-R-genome.xlsx"))
    gc()
    
    invisible(ca)
}

CalcNeoExpressedGenes <- function(ca) {
    exp.AA <- ca$EXP$genes$AA
    exp.RR <- ca$EXP$genes$RR
    exp.IA <- rownames(ca$fpkm$IA[apply(ca$fpkm$IA, 1, mean) > FPKM_CUTOFF, ])
    exp.IR <- rownames(ca$fpkm$IR[apply(ca$fpkm$IR, 1, mean) > FPKM_CUTOFF, ])
    exp.all <- ca$EXP$genes$ALL
    
    neo.exp.IA.from.AA <- setdiff(exp.IA, exp.AA)
    neo.exp.IA.from.RR <- intersect(neo.exp.IA.from.AA, exp.RR)
    neo.exp.IR.from.RR <- setdiff(exp.IR, exp.RR)
    neo.exp.IR.from.AA <- intersect(neo.exp.IR.from.RR, exp.AA)
    exp.IA.and.AA <- intersect(exp.AA, exp.IA)
    exp.IR.and.RR <- intersect(exp.RR, exp.IR)
    
    rep.exp.IA.from.AA <- setdiff(exp.AA, exp.IA)
    rep.exp.IR.from.RR <- setdiff(exp.RR, exp.IR)
    
    png(paste0(path.lib, "/ExpNeo-genes-venn.png"), 500, 500)
    v <- venn(list(AA = exp.AA, IA = exp.IA, IR = exp.IR, RR = exp.RR))
    dev.off()
    
    genes <- gofull <- gosimple <- vector("list", length = length(attr(v, "intersections")) + 1)
    names(genes) <- names(gofull) <- names(gosimple) <- c(names(attr(v, "intersections")), "II")
    for (vn in names(attr(v, "intersections"))) {
        .g <- attr(v, "intersections")[[vn]]
        e <- try(.go <- DoGO2(sig.genes = .g, all.genes = exp.all))
        if (class(e) == "try-error") {
            gofull[[vn]] <- data.frame(Desc = "no results.")
            gosimple[[vn]] <- data.frame(Desc = "no results.")
        } else {
            gofull[[vn]] <- .go$GOTABLE
            gosimple[[vn]] <- .go$GOSMPLTABLE
        }
        genes[[vn]] <- BindTairName(.g)
    }
    
    .go <- DoGO2(sig.genes = exp.ii.union, all.genes = exp.all)
    genes[["II"]] <- BindTairName(exp.ii.union)
    gofull[["II"]] <- .go$GOTABLE
    gosimple[["II"]] <- .go$GOSMPLTABLE
    
    names(gofull)    <- gsub(":", "-", names(gofull))
    names(gosimple) <- gsub(":", "-", names(gosimple))
    names(genes)     <- gsub(":", "-", names(genes))
    SaveExcel(gofull, paste0(path.lib, "/ExpNeo-genes-GO-full.xlsx"))
    gc()
    SaveExcel(gosimple, paste0(path.lib, "/ExpNeo-genes-GO-simple.xlsx"))
    gc()
    SaveExcel(genes, paste0(path.lib, "/ExpNeo-genes.xlsx"))
    gc() 
    
    invisible(ca)
}


IdentifyVEH <- function(ca) {
    .calchighcvgenes <- function(f, tag) {
        log10.f <- log10(f + 1)
        sd.log10.f <- apply(log10.f, 1, sd, na.rm = T)
        mu.log10.f <- apply(log10.f, 1, mean, na.rm = T)
        cv.log10.f <- sd.log10.f / mu.log10.f
        
        g <- ggplot(data.frame(cv = cv.log10.f, mean = mu.log10.f, gene = ifelse( (cv.log10.f > 0.20) & (mu.log10.f > 1.0), 'VEG', 'nonVEG')),
                    aes(x = mean, y = cv, color = gene))
        g <- g + geom_point(alpha = 0.2, size = 2) + scale_color_manual(values =  c('#999999', '#E41A1C'))
        g <- g + xlab('mean(log10FPKM)') + ylab('cv(log10FPKM)')
        png(paste0(path.expbias, '/CV.disp.', tag, '.png'), 500, 400)
        plot(g)
        dev.off()
        
        highcv <- (cv.log10.f > 0.20) & (mu.log10.f > 1.0)
        highcv.genes <- names(highcv)[highcv]
        goobj <- DoGO2(sig.genes = highcv.genes, all.genes = ca$EXP$genes$ALL)
        SaveExcel(list(full=goobj$GOTABLE, simlify=goobj$GOSMPLTABLE),
                  file.name = paste0(path.expbias, "/GO_", tag, "_CV0.20_and_mu1.0.xls"))
        invisible(list(go = goobj$GOSMPLTABLE, gene = highcv.genes,
                       df = data.frame(gene = highcv.genes, sd = sd.log10.f[highcv.genes],
                                  mean = mu.log10.f[highcv.genes], cv = cv.log10.f[highcv.genes])))
    }
    
    A <- ca$fpkm$AA[ca$EXP$genes$ALL, ]
    R <- ca$fpkm$RR[ca$EXP$genes$ALL, ]
    I <- ca$fpkm$II[ca$EXP$genes$ALL, ]
    a <- ca$fpkm$IA[ca$EXP$genes$ALL, ]
    r <- ca$fpkm$IR[ca$EXP$genes$ALL, ]
    A.go <- .calchighcvgenes(A, "AA")
    R.go <- .calchighcvgenes(R, "RR")
    a.go <- .calchighcvgenes(a, "IA")
    r.go <- .calchighcvgenes(r, "IR")
    I.go <- .calchighcvgenes(I, "II")
    print(paste0("# highCV genes in A: ", length(A.go$gene)))
    print(paste0("# highCV genes in R: ", length(R.go$gene)))
    print(paste0("# highCV genes in a: ", length(a.go$gene)))
    print(paste0("# highCV genes in r: ", length(r.go$gene)))
    print(paste0("# highCV genes in I: ", length(I.go$gene)))
    
    vlist <- list(IA = a.go$gene, IR = r.go$gene)#), A = A.go$gene, R = R.go$gene)
    
    highcv.genes.df <- list(Sheet1 = BindTairName(a.go$df),
                            Sheet2 = BindTairName(r.go$df),
                            Sheet3 = BindTairName(A.go$df),
                            Sheet4 = BindTairName(R.go$df))
    SaveExcel(highcv.genes.df, file.name = paste0(path.expbias, "/VEHs.IA_IR_AA_RR.xlsx"))
 
    goterms <- unique(c(a.go$go$ID[a.go$go$qvalue < 0.1],
                        r.go$go$ID[r.go$go$qvalue < 0.1],
                        A.go$go$ID[A.go$go$qvalue < 0.1],
                        R.go$go$ID[R.go$go$qvalue < 0.1]))
    go2desc <- rep("", length = length(goterms))
    gomat <- matrix(NA, ncol = 4, nrow = length(goterms))

    rownames(gomat) <- names(go2desc) <- goterms
    colnames(gomat) <- c("IA", "IR", "A", "R")
    go2desc[a.go$go$ID[a.go$go$qvalue < 0.1]] <- a.go$go$Description[a.go$go$qvalue < 0.1]
    go2desc[r.go$go$ID[r.go$go$qvalue < 0.1]] <- r.go$go$Description[r.go$go$qvalue < 0.1]
    go2desc[A.go$go$ID[A.go$go$qvalue < 0.1]] <- A.go$go$Description[A.go$go$qvalue < 0.1]
    go2desc[R.go$go$ID[R.go$go$qvalue < 0.1]] <- R.go$go$Description[R.go$go$qvalue < 0.1]

    gomat[a.go$go$ID[a.go$go$qvalue < 0.1], "IA"] <- -log10(a.go$go$qvalue[a.go$go$qvalue < 0.1])
    gomat[r.go$go$ID[r.go$go$qvalue < 0.1], "IR"] <- -log10(r.go$go$qvalue[r.go$go$qvalue < 0.1])
    gomat[A.go$go$ID[A.go$go$qvalue < 0.1], "A"] <- -log10(A.go$go$qvalue[A.go$go$qvalue < 0.1])
    gomat[R.go$go$ID[R.go$go$qvalue < 0.1], "R"] <- -log10(R.go$go$qvalue[R.go$go$qvalue < 0.1])

    df <- data.frame(gomat, desc = go2desc)
    SaveExcel(list(Sheet1 = df), file.name = paste0(path.expbias, "/VEH_GOresults.xlsx"))

    
    ca
}


plotFCRsquared <- function(ca) {
    A <- log10(ca$fpkm$AA + 1)[ca$EXP$genes$ALL, ]
    R <- log10(ca$fpkm$RR + 1)[ca$EXP$genes$ALL, ]
    a <- log10(ca$fpkm$IA + 1)[ca$EXP$genes$ALL, ]
    r <- log10(ca$fpkm$IR + 1)[ca$EXP$genes$ALL, ]
    A.fc <- (A - A[, 1])[, -1]
    R.fc <- (R - R[, 1])[, -1]
    a.fc <- (a - a[, 1])[, -1]
    r.fc <- (r - r[, 1])[, -1]
    
    f <- log10(ca$fpkm$II + 1)[ca$EXP$genes$ALL, ]
    sd.log10.f <- apply(f, 1, sd, na.rm = T)
    mu.log10.f <- apply(f, 1, mean, na.rm = T)
    cv.log10.f <- sd.log10.f / mu.log10.f
    
    fcmat <- cbind(A.fc, a.fc, r.fc, R.fc)
    fcmat <- fcmat[cv.log10.f > 0.6, ]
    colnames(fcmat) <- paste0(rep(c("A", "IA", "IR", "R"), each = 8),
                              " (", TIMELABELS[-1], ")")
    
    rsqmatrix <- matrix(NA, ncol = ncol(fcmat), nrow = ncol(fcmat))
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- colnames(fcmat)
    for (ci in 1:ncol(rsqmatrix)) {
        for (ri in 1:nrow(rsqmatrix)) {
            rsqmatrix[ri, ci] <- summary(lm(fcmat[, ri] ~ fcmat[, ci]))$r.squared
        }
    }
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- rep(TIMELABELS[-1], times = 4)

    png(paste0(path.de, "/log10FC.fpkm.Rsquared.png"), 800, 800)
    heatmap.2(rsqmatrix, trace = "none", scale="none", dendrogram="none", Colv = F, Rowv = F,
              density.info = "none", col = COLSFUNC(), margins = c(4, 4),
              key.xlab = "R-squared value", key.title = "", sepcolor = "white",
              colsep = c(1:4 * 8), rowsep = c(1:4 * 8), sepwidth = c(0.1, 0.1))
    dev.off()
}



plotScatter <- function(ca, g, path.prefix) {
    for (.g in g) {
        .t <- CARHR2TAIR[[.g]]
        .s <- CARHR2NAME[[.g]]
        if (is.null(.t)) {
            gname <- .g
        } else {
            if (is.null(.s)) {.s <- ""
            } else { .s <- paste0(" / ", .s)}
            gname <- paste0(.g, " (", .t, .s, ")")
            
        }
        exp.A <- log10(ca$fpkm$AA[.g, ] + 1)
        exp.a <- log10(ca$fpkm$IA[.g, ] + 1)
        exp.r <- log10(ca$fpkm$IR[.g, ] + 1)
        exp.R <- log10(ca$fpkm$RR[.g, ] + 1)
        exp.I <- log10(ca$fpkm$IA[.g, ] * (1/2) + ca$fpkm$IA[.g, ] * (2/2) + 1)
        ARatio <- ca$aorigin.ratio[.g, ]
        lb <- gsub(" hr", "", TIMELABELS)
        png(paste0(path.prefix, "/", .g, ".png"), 730, 320)
        plot(0, 0, type = "n",  xlim = c(-0.2, 1.2), axes = F,
             ylim = c(-0.2, max(3, max(c(exp.a, exp.r, exp.A, exp.R)))) + 0.2,
             xlab = "A-origin ratio", ylab = "log10FPKM",
             main = gname)
        grid()
        axis(1, c(0, 0.2, 1/3, 0.4, 0.6, 0.8, 1.0), c("0", "0.2", "1/3", "0.4", "0.6", "0.8", "1.0"))
        axis(2)
        #abline(v = c(0, 1/3, 2/3, 1), col = "#6f6f6f", lty = 2)
        abline(v = c(0, 1), col = "#000000", lty = 1)
        text(rev(-0.02 * 1:9), exp.R, labels = lb, col = brewer.pal(8, "Dark2"))
        lines(rev(-0.02 * 1:9), exp.R, col = "#999999")
        text(1 + 0.02*1:9, exp.A, labels = lb, col = brewer.pal(8, "Dark2"))
        lines(1 + 0.02 * 1:9, exp.A, col = "#999999")
        
        lines(ARatio, exp.I, col = "#999999", lwd = 1.4)
        text(ARatio, exp.I, labels = lb, col = brewer.pal(8, "Dark2"), cex = 2)
        dev.off()
    }
    

}






DoGO <- function(sig.genes = NULL, all.genes = NULL, p.values = NULL, ontology = "BP",
                 p.cutoff = 1, q.cutoff = 0.1, method = "over-representation") {
    all.genes <- as.character(all.genes)
    sig.genes <- as.character(sig.genes)
    all.genes.at <- as.character(unlist(CARHR2TAIR[all.genes]))
    sig.genes.at <- as.character(unlist(CARHR2TAIR[sig.genes]))
    all.genes.at <- all.genes.at[all.genes.at != ""]
    sig.genes.at <- sig.genes.at[sig.genes.at != ""]
    ego <- enrichGO(gene = sig.genes.at, universe = all.genes.at, OrgDb = org.At.tair.db,
                    ont = ontology, pAdjustMethod = "BH", keyType = "TAIR",
                    pvalueCutoff = p.cutoff, qvalueCutoff = q.cutoff, readable = FALSE)
    rgo <- summary(ego)
    list(TABLE = rgo, OBJ = ego)
}




BindTairName <- function(dat) {
    if (!is.matrix(dat) && !is.data.frame(dat)) {
        dat <- as.matrix(dat)
        colnames(dat) <- "CARHR"
        rownames(dat) <- dat
    }
    rnm <- colnames(dat)
    tairid <- as.character((CARHR2TAIR[as.character(rownames(dat))]))
    genename <- as.character((CARHR2NAME[as.character(rownames(dat))]))
    genedesc <- as.character((CARHR2DESC[as.character(rownames(dat))]))
    tairid[tairid == "NULL"] <- ""
    genename[genename == "NULL"] <- ""
    genedesc[genedesc == "NULL"] <- ""
    dat <- data.frame(dat, TAIR = tairid, name = genename, description = genedesc, stringsAsFactors = F)
    dat[is.na(dat)] <- ""
    colnames(dat) <- c(rnm, "TAIR", "name", "description")
    dat
}








PlotRsqHeatmap <- function(f, type = "all") {
    A <- log10(ca$fpkm$AA + 1)
    R <- log10(ca$fpkm$RR + 1)
    a <- log10(ca$fpkm$IA + 1)
    r <- log10(ca$fpkm$IR + 1)
    y <- cbind(A, a, r, R)
    
    I <- log10(ca$fpkm$II + 1)
    cv.I <- apply(I, 1, sd) / apply(I, 1, mean)
    
    y <- y[cv.I > 0.8, ]
    rsqmatrix <- NULL
    rsqmatrix <- matrix(NA, ncol = 4 * 9, nrow = 4 * 9)
    colnames(rsqmatrix) <- rownames(rsqmatrix) <- paste0(rep(SPLABELS, each = 9), " (", rep(TIMELABELS, times = 4), ")")
    for (ci in 1:ncol(rsqmatrix)) {
        for (ri in 1:nrow(rsqmatrix)) {
            rsqmatrix[ri, ci] <- summary(lm(y[, ri] ~ y[, ci]))$r.squared
        }
    }
    heatmap.2(rsqmatrix, trace = "none", scale="none", dendrogram="none", Colv = F, Rowv = F,
              density.info = "none", col = COLSFUNCC2(), margins = c(10, 10),
              key.xlab = "R-squared value", key.title = "", sepcolor = "white",
              colsep = c(1:4 * 9), rowsep = c(1:4 * 9), sepwidth = c(0.1, 0.1))
}





#SimplifyGO <- function(gotable) {
#    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
#    gotable <- gotable[intersect(gotable$ID, target.terms), ]
#    gotable$qvalue <- p.adjust(gotable$pvalue, method = "BH")
#    gotable
#}















.common.RData  <- paste0(PATH$lib, "/src/common.RData")
.ca.RData <- paste0(PATH$lib, "/src/ca.RData")
load(.common.RData)
load(.ca.RData)
stop()


if (file.exists(.ca.RData)) {
    load(.ca.RData)
} else {
    design <- GetDesign()
    ca <- NULL
    ca <- GetCounts(ca, design)
    ca <- CalcFPKM(ca)
    ca <- CalcExpressedGenes(ca)
    #ca <- CalcNeoExpressedGenes(ca)
    ca <- IdentifyVEH(ca)
    ca <- PlotAoriginRatioDist(ca)
    save(ca, design, file = .ca.RData)
}




## scatterplot of A-origin ratio between 0 hr and 8 hr
fa <- ca$fpkm$IA[, c(1, 4)]
fr <- ca$fpkm$IR[, c(1, 4)]
f <- cbind(fa, fr)
f <- rownames(f[apply(f, 1, mean) > 1, ])
ardf <- data.frame(H00 = ca$aorigin.ratio[f, 1],
                   H08 = ca$aorigin.ratio[f, 4]) # 0,2,4,8
ardf <- ardf[rowSums(is.na(ardf)) == 0, ]
g <- ggplot(ardf, aes(H00, H08))
g <- g + geom_point(size = 1, alpha = 0.5)
g <- g + theme_bw() + xlab("0 hr") + ylab("8 hr")
pdf(paste0(path.expbias, "/0h8h-aorigin-ratio.scatter.pdf"), 4.5, 4.5)
ggMarginal(g, type = "histogram", bins = 40, col = "#363636", fill ="#363636")
dev.off()


## number of newly-expressed genes and silenced genes
lapply(ca$EXP$genes, length)
# A - Ia
length(intersect(ca$EXP$genes$AA, ca$EXP$genes$IA)) / length(union(ca$EXP$genes$AA, ca$EXP$genes$IA))
# R - Ir
length(intersect(ca$EXP$genes$RR, ca$EXP$genes$IR)) / length(union(ca$EXP$genes$RR, ca$EXP$genes$IR))
# A - R
length(intersect(ca$EXP$genes$AA, ca$EXP$genes$RR)) / length(union(ca$EXP$genes$AA, ca$EXP$genes$RR))
# Ia - Ir
length(intersect(ca$EXP$genes$IA, ca$EXP$genes$IR)) / length(union(ca$EXP$genes$IA, ca$EXP$genes$IR))




## gene expression profiles (z-scored)
g <- c('CARHR049560', # STM/BUM
       'CARHR279690', # BAM1
       'CARHR158740', # BAM2
       'CARHR227100', # BAM3
       'CARHR202780', # PLT3
       'CARHR274370') # REV
g.fpkm <- cbind(ca$fpkm$AA[g, ], ca$fpkm$IA[g, ], ca$fpkm$IR[g, ], ca$fpkm$RR[g, ])
rownames(g.fpkm) <- c('STM', 'BAM1', 'BAM2', 'BAM3', 'PLT3', 'REV')
colnames(g.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                           '  ', colnames(g.fpkm))
.g.fpkm <- t(apply(log10(g.fpkm + 1), 1, scale))
colnames(.g.fpkm) <- colnames(g.fpkm)
g.fpkm <- .g.fpkm
g.fpkm.df <- melt(g.fpkm)
colnames(g.fpkm.df) <- c('gene', 'library', 'zscore')
ghm <- ggplot(g.fpkm.df, aes(x = library, y = gene, fill = zscore))
ghm <- ghm + geom_tile() + scale_fill_gradientn('value',
            colours = rev(brewer.pal(9, 'Spectral')))
ghm <- ghm + xlab('') + ylab('') + coord_fixed()
ghm <- ghm + theme_bw() + guides(fill = guide_legend(title = 'z-score'))
ghm <- ghm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   strip.background = element_rect(fill = "white", colour = "white"))
pdf(paste0(path.geneprof, '/meristem.selected.homeologs.pdf'), 10.8, 3)
plot(ghm)
dev.off()






## meristem associated genes (gene expression)
meristem.genes <- unique(c(GO2CARHR[['GO:0035266']], GO2CARHR[['GO:0010075']], GO2CARHR[['GO:0048509']]))
meristem.genes <- intersect(meristem.genes, rownames(ca$fpkm$AA))
meristem.fpkm <- cbind(ca$fpkm$AA[meristem.genes, ], ca$fpkm$IA[meristem.genes, ], 
                      ca$fpkm$IR[meristem.genes, ], ca$fpkm$RR[meristem.genes, ])
colnames(meristem.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                           '  ', colnames(meristem.fpkm))
meristem.fpkm <- log10(meristem.fpkm + 1)
meristem.fpkm <- BindTairName(meristem.fpkm)


## water deprivation
watdep.genes <- GO2CARHR[['GO:0009414']]
watdep.genes <- intersect(watdep.genes, rownames(ca$fpkm$AA))
watdep.fpkm <- cbind(ca$fpkm$AA[watdep.genes, ], ca$fpkm$IA[watdep.genes, ], 
                     ca$fpkm$IR[watdep.genes, ], ca$fpkm$RR[watdep.genes, ])
colnames(watdep.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(watdep.fpkm))
watdep.fpkm <- log10(watdep.fpkm + 1)
watdep.fpkm <- BindTairName(watdep.fpkm)


## response to ethylene
ethylene.genes <- GO2CARHR[['GO:0009723']]
ethylene.genes <- intersect(ethylene.genes, rownames(ca$fpkm$AA))
ethylene.fpkm <- cbind(ca$fpkm$AA[ethylene.genes, ], ca$fpkm$IA[ethylene.genes, ], 
                       ca$fpkm$IR[ethylene.genes, ], ca$fpkm$RR[ethylene.genes, ])
colnames(ethylene.fpkm) <- paste0(rep(c('A', 'IA', 'IR', 'R'), each = 9),
                                '  ', colnames(ethylene.fpkm))
ethylene.fpkm <- log10(ethylene.fpkm + 1)
ethylene.fpkm <- BindTairName(ethylene.fpkm)

SaveExcel(list(meristem = meristem.fpkm, waterdep = watdep.fpkm, ethylene = ethylene.fpkm),
          file = paste0(path.geneprof, '/gene.prof.xls'))









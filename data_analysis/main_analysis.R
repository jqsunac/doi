##
##
##

set.seed(2016)
suppressMessages(require(MASS))
suppressMessages(require(RUnit))
suppressMessages(require(knitr))
suppressMessages(require(XLConnect))
suppressMessages(require(RColorBrewer))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))
suppressMessages(require(ggExtra))
suppressMessages(require(gplots))
suppressMessages(require(VennDiagram))
suppressMessages(require(grid))
suppressMessages(require(scales))
suppressMessages(require(outliers))
suppressMessages(require(fBasics))
suppressMessages(require(edgeR))
suppressMessages(require(TCC))
suppressMessages(require(GO.db))
suppressMessages(require(org.At.tair.db))
suppressMessages(require(clusterProfiler))
suppressMessages(require(plyr))
suppressMessages(require(genefilter))
suppressMessages(require(googleVis))
suppressMessages(require(venneuler))
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
    flx  = .brewer.set1[3],     # C. flexuosa
    flxA = .brewer.set1[4],     # C. flexuosa (A)
    flxH = .brewer.set1[5],     # C. flexuosa (H)
    DE   = .brewer.set1[5],     # Differentially expressed genes
    NDE  = .brewer.set1[9]      # Non-differentially expressed genes
)

COLSFUNC <- function(n = 100) {
    #.colfunc <- colorRampPalette(c("#FFFFFF", "#D6604D", "#67000D"))
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu")))
    .colfunc(n)
}
COLSFUNC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(c("#053061", "#4393C3", "#FFFFFF", "#D6604D", "#67001F"))
    .colfunc(n)
}
COLSFUNCC1 <- function(n = 100) {
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu")))
    .colfunc(n)
}
COLSFUNCC2 <- function(n = 100) {
    .colfunc <- colorRampPalette(rev(brewer.pal(9,"RdBu"))[5:9])
    .colfunc(n)
}
COLSFUNC23 <- function(n = 100) {
    .colfunc <- colorRampPalette(c("#053061", "#FFFFFF", "#D6604D", "#67001F"))
    .colfunc(n)
}




SPLABELS  <- c("AA", "IA", "IR", "RR")
SPLABELS2 <- c("AA", "IA", "II", "IR", "RR")

TIMELABELS  <- c("0 hr", "2 hr", "4 hr", "8 hr", "12 hr", "24 hr", "48 hr", "72 hr", "96 hr")
TIMELABELS2 <-  c("0 hr vs. 2 hr",   "0 hr vs. 4 hr",  "0 hr vs. 8 hr", "0 hr vs. 12 hr",
                  "0 hr vs. 24 hr", "0 hr vs. 48 hr", "0 hr vs. 72 hr", "0 hr vs. 96 hr")
INCH2CM <- 1 / 2.54


HASH6 <- list(`pat` = "CaA", `mat` = "CrR", `hyb` = "CiI", `hybP` = "CiA", `hybM` = "CiR")
HASH7 <- list(`ama` = "CaA", `riv` = "CrR", `ins` = "CiI", `insA` = "CiA", `insR` = "CiR")
HASH2 <- list(`1000` = "CaA",   `0100` = "CiA",   `0010` = "CiR",    `0001` = "CrR",
              `1101` = "CaA:CiA:CrR", `1011` = "CaA:CiR:CrR",
              `1110` = "CaA:CiA:CiR", `0111` = "CiA:CiR:CrR",
              `1111` = "CaA:CiA:CiR:CrR", `1001` = "CaA:CrR",
              `1100` = "CaA:CiA",  `1010` = "CaA:CiR",
              `0101` = "CiA:CrR",  `0011` = "CiR:CrR", `0110` = "CiA:CiR")
HASH4 <- list(`CaA` = "1000",  `CiA` =  "0100",  `CiR` = "0010",  `CrR` = "0001",
              `CaA:CiA:CrR` = "1101", `CaA:CiR:CrR` = "1011",
              `CaA:CiA:CiR` = "1110", `CiA:CiR:CrR` = "0111",
              `CaA:CiA:CiR:CrR` = "1111", `CaA:CrR` = "1001",
              `CaA:CiA` = "1100", `CaA:CiR` = "1010",
              `CiA:CrR` = "0101", `CiR:CrR` = "0011", `CiA:CiR` = "0110")





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
    workbook <- loadWorkbook(file.name, create = TRUE)
    createSheet(workbook, names(dat))
    writeWorksheet(workbook, dat, names(dat), header = TRUE)
    saveWorkbook(workbook)
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
    gene.names  <- gene.names[-c(grep("^__", gene.names))]
    
    counts.AA <- counts.aorigin[, 1:9] + counts.rorigin[, 1:9] + counts.common[, 1:9]
    counts.RR <- counts.aorigin[, 19:27] + counts.rorigin[, 19:27] + counts.common[, 19:27]
    
    aorigin.ratio <- (counts.aorigin / (counts.aorigin + counts.rorigin))[, 10:18]
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
    .counts.I <- ca$counts$IA + ca$counts$IR
    cpm.I <- sweep(.counts.I, 2, 1e6 / colSums(.counts.I), "*")
    .fpkm.I <- sweep(cpm.I, 1, 1e3 / gene.length[, 1], "*")
    fpkm$IA <- .fpkm.I * ca$aorigin.ratio
    fpkm$IR <- .fpkm.I * (1 - ca$aorigin.ratio)
    fpkm$IA[is.na(fpkm$IA)] <- 0
    fpkm$IR[is.na(fpkm$IR)] <- 0
    
    ca$gene.length <- gene.length
    ca$fpkm <- fpkm
    ca
}



RunEdgeRPairewise <- function(c1, c2, f1, f2) {
    keep <- (f1 > FPKM_CUTOFF | f2 > FPKM_CUTOFF)
    cnt <- as.matrix(round(cbind(c1, c2)))[keep, ]
    cnt  <- cnt[(rowSums(cnt) > 0), ]
    t <- new("TCC", cnt, group = c("A", "B"))
    t <- calcNormFactors(t)
    d <- DGEList(counts = cnt, group = c("A", "B"), norm.factors = t$norm.factors)
    d <- estimateGLMCommonDisp(d, method = "deviance", robust = T, subset = NULL)
    res <- exactTest(d)
    res.df <- as.data.frame(topTags(res, n = nrow(cnt)))
    res.df <- data.frame(gene = as.character(rownames(res.df)), res.df[, c(1, 2, 3, 4)])
    colnames(res.df) <- c("ID", "logFC", "logCPM", "p.value", "q.value")
    res.df
}



IdentifyDEGs <- function(ca, p.cutoff = NULL, q.cutoff = NULL) {
    if (is.null(p.cutoff) && is.null(q.cutoff)) stop("One of 'p.cutoff' and 'q.cutoff' should be given.")
    countslist <- ca$counts
    fpkmlist   <- ca$fpkm
    fpkmlist$IA <- fpkmlist$IR <- ca$fpkm$IA + ca$fpkm$IR
    
    tmpl.list <- vector("list", length(ca$counts))
    names(tmpl.list) <- names(ca$counts)

    ca$DEG$LIST <- ca$DEG$UP <- ca$DEG$FL <- ca$DEG$DW <- tmpl.list
    ca$DEG$SUM <- ca$DEG$SUMP <- matrix(0, nrow = 2 * length(ca$DEG$LIST), ncol = ncol(ca$counts[[1]]) - 1)
    colnames(ca$DEG$SUM) <- colnames(ca$DEG$SUMP) <- TIMELABELS[-1]
    rownames(ca$DEG$SUM) <- rownames(ca$DEG$SUMP) <- paste(rep(names(ca$DEG$LIST), each = 2), c("up","dw"), sep = ".")
    
    ma.coordinates <- vector("list", length(ca$DEG$LIST))
    names(ma.coordinates) <- names(ca$DEG$LIST)

    for (i in 1:length(ca$DEG$LIST)) {
        ma.coordinates[[i]] <- ca$DEG$LIST[[i]] <-
            ca$DEG$UP[[i]] <- ca$DEG$FL[[i]] <- ca$DEG$DW[[i]] <-
            vector("list", ncol(countslist[[i]]) - 1)
        names(ma.coordinates[[i]]) <- names(ca$DEG$LIST[[i]]) <-
            names(ca$DEG$UP[[i]]) <- names(ca$DEG$FL[[i]]) <- names(ca$DEG$DW[[i]]) <-
            colnames(ca$DEG$SUM)
        for (j in 1:(ncol(ca$counts[[i]]) - 1)) {
            message("[IdentifyDEGs] Identying DEGs between libraries of ",
                    colnames(countslist[[i]])[1], " and ", colnames(countslist[[i]])[j + 1], ".")

            ca$DEG$LIST[[i]][[j]] <- RunEdgeRPairewise(countslist[[i]][, 1], countslist[[i]][, j + 1],
                                                       fpkmlist[[i]][, 1],   fpkmlist[[i]][, j + 1])
            if (!is.null(p.cutoff)) keeped.ids <- (ca$DEG$LIST[[i]][[j]]$p.value < p.cutoff)
            if (!is.null(q.cutoff)) keeped.ids <- (ca$DEG$LIST[[i]][[j]]$q.value < q.cutoff)
            ca$DEG$UP[[i]][[j]] <- as.character(ca$DEG$LIST[[i]][[j]]$ID[keeped.ids & (ca$DEG$LIST[[i]][[j]]$logFC > 0)])
            ca$DEG$DW[[i]][[j]] <- as.character(ca$DEG$LIST[[i]][[j]]$ID[keeped.ids & (ca$DEG$LIST[[i]][[j]]$logFC < 0)])
            ca$DEG$FL[[i]][[j]] <- as.character(setdiff(as.character(ca$DEG$LIST[[i]][[j]]$ID),
                                                union(ca$DEG$UP[[i]][[j]], ca$DEG$DW[[i]][[j]])))
            ca$DEG$SUM[(i - 1) * 2 + 1, j] <- length(ca$DEG$UP[[i]][[j]])
            ca$DEG$SUM[(i - 1) * 2 + 2, j] <- length(ca$DEG$DW[[i]][[j]])
            ca$DEG$SUMP[(i - 1) * 2 + 1, j] <- length(ca$DEG$UP[[i]][[j]]) / nrow(ca$DEG$LIST[[i]][[j]]) * 100
            ca$DEG$SUMP[(i - 1) * 2 + 2, j] <- length(ca$DEG$DW[[i]][[j]]) / nrow(ca$DEG$LIST[[i]][[j]]) * 100
            ma.coordinates[[i]][[j]] <- data.frame(gene = ca$DEG$LIST[[i]][[j]]$ID, Type = ifelse(keeped.ids, "DEG", "NDEG"),
                                                   logFC = ca$DEG$LIST[[i]][[j]]$logFC, logCPM = ca$DEG$LIST[[i]][[j]]$logCPM)

            message("[IdentifyDEGs] Identified up-regulated DEGs: ", length(ca$DEG$UP[[i]][[j]]), " .")
            message("[IdentifyDEGs] Identified down-regulated DEGs: ", length(ca$DEG$DW[[i]][[j]]), " .")
        }
        label_1 <- colnames(ca$counts[[i]])
        names(ca$DEG$LIST[[i]]) <- paste("00", gsub("[a-zA-Z ]", "", label_1[-1]), sep = "_")
        names(ca$DEG$UP[[i]]) <- names(ca$DEG$FL[[i]]) <- names(ca$DEG$DW[[i]]) <- gsub("[a-zA-Z.]", "", names(ca$DEG$LIST[[i]]))
    }
    gdatdf <- NULL
    for (i in 1:length(ma.coordinates)) {
        for (j in 1:length(ma.coordinates[[i]])) {
            gdatdf <- rbind(gdatdf,
                            data.frame(Species = names(ma.coordinates)[i],
                                       Time    = names(ma.coordinates[[i]])[j],
                                       ma.coordinates[[i]][[j]]))
        }
    }

    gdatdf$Species <- as.character(gdatdf$Species)

    gdatdf$Species <- factor(gdatdf$Species, levels = SPLABELS2)
    gdatdf$Time    <- factor(gdatdf$Time,    levels = TIMELABELS[-1])
    gg <- ggplot(gdatdf, aes(x = logCPM, y = logFC, color = Type))
    gg <- gg + geom_point(size = 0.1)
    gg <- gg + theme_bw()
    gg <- gg + theme(legend.position = "none", strip.background = element_rect(fill = "white", colour = "white"))
    gg <- gg + scale_colour_manual(values = c(COLS$DE, COLS$NDE))
    gg <- gg + facet_grid(Species ~ Time)
    gg <- gg + ylab("logFC") + xlab("logCPM")
    png(paste0(path.de, "/maplots.png"), 700, 600)
    plot(gg)
    dev.off()
    
    ca$DEG$plots$maplots <- gg
    ca$DEG$stats$deg.sump <- ca$DEG$SUMP
    ca$DEG$stats$deg.sum  <- ca$DEG$SUM
    ca$DEG$used.cutoff <- list(FPKM = FPKM_CUTOFF, p = p.cutoff, q = q.cutoff)
    ca
}







IdentifyDEGsFCPlot <- function(ca, p.cutoff = NULL, q.cutoff = NULL) {
    message("Starting to identify DEGs.")

    if (is.null(p.cutoff) && is.null(q.cutoff)) stop("One of 'p.cutoff' and 'q.cutoff' should be given.")

    tmpl.list <- vector("list", length(ca$parsed.counts))
    names(tmpl.list) <- names(ca$parsed.counts)

    degsheet <- NULL

    for (sp in names(ca$DEG$LIST)) {
        for (tm in names(ca$DEG$LIST[[sp]])) {
            is.deg <- as.logical(ca$DEG$LIST[[sp]][[tm]]$q.value < q.cutoff)
            up.fc0 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC > 0)
            up.fc4 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC > 2)
            up.fc8 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC > 3)
            dw.fc0 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC < -0)
            dw.fc4 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC < -2)
            dw.fc8 <- as.logical(ca$DEG$LIST[[sp]][[tm]]$logFC < -3)

            degsheet.1 <- data.frame(class = sp, time = tm, deg = "up-regulated", foldchange = "2 < FC ≤ 4",
                                     value = sum(is.deg & up.fc0) - sum(is.deg & up.fc4))
            degsheet.2 <- data.frame(class = sp, time = tm, deg = "up-regulated", foldchange = "4 < FC ≤ 8",
                                     value = sum(is.deg & up.fc4) - sum(is.deg & up.fc8))
            degsheet.3 <- data.frame(class = sp, time = tm, deg = "up-regulated", foldchange = "8 < FC",
                                     value = sum(is.deg & up.fc8))
            degsheet.4 <- data.frame(class = sp, time = tm, deg = "down-regulated", foldchange = "2 < FC ≤ 4",
                                     value = sum(is.deg & dw.fc0) - sum(is.deg & dw.fc4))
            degsheet.5 <- data.frame(class = sp, time = tm, deg = "down-regulated", foldchange = "4 < FC ≤ 8",
                                     value = sum(is.deg & dw.fc4) - sum(is.deg & dw.fc8))
            degsheet.6 <- data.frame(class = sp, time = tm, deg = "down-regulated", foldchange = "8 < FC",
                                     value = sum(is.deg & dw.fc8))
            degsheet <- rbind(degsheet, degsheet.1, degsheet.2, degsheet.3,
                                        degsheet.4, degsheet.5, degsheet.6)
        }
    }

    degsheet$class <- factor(degsheet$class, levels = SPLABELS)
    degsheet$time <- as.character(degsheet$time)
    degsheet$time <- gsub("00_12", TIMELABELS[5], degsheet$time)
    degsheet$time <- gsub("00_24", TIMELABELS[6], degsheet$time)
    degsheet$time <- gsub("00_48", TIMELABELS[7], degsheet$time)
    degsheet$time <- gsub("00_72", TIMELABELS[8], degsheet$time)
    degsheet$time <- gsub("00_96", TIMELABELS[9], degsheet$time)
    degsheet$time <- gsub("00_2", TIMELABELS[2], degsheet$time)
    degsheet$time <- gsub("00_4", TIMELABELS[3], degsheet$time)
    degsheet$time <- gsub("00_8", TIMELABELS[4], degsheet$time)
    degsheet$time <- factor(degsheet$time, levels = TIMELABELS[-1])
    degsheet$deg <- factor(degsheet$deg, levels = c("up-regulated", "down-regulated"))
    degsheet$foldchange <- factor(degsheet$foldchange, levels = c("2 < FC ≤ 4", "4 < FC ≤ 8", "8 < FC"))

    gg <- ggplot(degsheet, aes(x = time, y = value, fill = class, alpha = foldchange))
    gg <- gg + geom_bar(stat = "identity")
    gg <- gg + theme_bw()
    gg <- gg + facet_grid(deg ~ class)
    gg <- gg + scale_fill_manual(values = c(COLS$ama, COLS$insA, COLS$insR, COLS$riv))
    gg <- gg + scale_alpha_manual(values=c(0.40, 0.70, 1))
    gg <- gg +  theme(axis.text.x = element_text(angle = 90, hjust = 1),
                      strip.background = element_rect(fill = "white", colour = "white"))
    gg <- gg + ylab("Time after submergence") + xlab("Number of DEGs")
    ca$DEG$plots$num.of.degs <- gg
    
    png(paste0(path.de, "/DEG_numberofdegs.png"), 500, 300)
    plot(gg)
    dev.off()
    
    ca
}



PlotFPKMScatter.plot_xy <- function(A = NULL, R = NULL, a = NULL, r = NULL, gid, gname, gdesc) {
    if (is.null(A) && is.null(R)) {
        A = ca$fpkm$AA[gid, ] * (1/3)
        R = ca$fpkm$RR[gid, ] * (2/3)
        a = ca$fpkm$IA[gid, ]
        r = ca$fpkm$IR[gid, ]
        gname <- ifelse(!is.null(CARHR2NAME[[gid]]), CARHR2NAME[[gid]], "")
        gdesc <- ifelse(!is.null(CARHR2DESC[[gid]]), CARHR2DESC[[gid]], "")
    }
    log10A <- log10(A + 1)
    log10a <- log10(a + 1)
    log10R <- log10(R + 1)
    log10r <- log10(r + 1)
    df <- rbind(data.frame(A = log10A, R = log10R, time = gsub(" hr", "", TIMELABELS), allopolyploid = "synthetic"),
                data.frame(A = log10a, R = log10r, time = gsub(" hr", "", TIMELABELS), allopolyploid = "C. insueta"))
    g <- ggplot(df, aes(x = A, y = R, label = time, colour = allopolyploid, group = allopolyploid))
    g <- g + geom_abline(intercept = 0, slope = 1, colour = "#AAAAAA", linetype = 2)
    g <- g + geom_abline(intercept = 0, slope = log10(3)/log10(2), colour = "#AAAAAA", linetype = 2)
    g <- g + geom_path() + theme_bw() + coord_fixed()
    g <- g + ggtitle(bquote(atop(.(paste0(gid, " ", gname)), atop(italic(.(gdesc)), "")))) 
    g <- g + geom_text(size = 4, colour = c(rep("#E41A1C", 9), rep("#377E88", 9))) 
    g <- g + xlab("log10(FPKM + 1) on A-genome") + ylab("log10(FPKM + 1) on R-genome")
    g <- g + scale_color_manual(values = c("#B3CDE3", "#FBB4AE"))
    g <- g + xlim(0, max(c(3, df$A, df$R))) + ylim(0, max(c(3, df$A, df$R)))
    g
}
    
PlotFPKMScatter.plot_rh <- function(A = NULL, R = NULL, a = NULL, r = NULL, gid, gname, gdesc) {
    if (is.null(A) && is.null(R)) {
        A = ca$fpkm$AA[gid, ] * (1/3)
        R = ca$fpkm$RR[gid, ] * (2/3)
        a = ca$fpkm$IA[gid, ]
        r = ca$fpkm$IR[gid, ]
        gname <- ifelse(!is.null(CARHR2NAME[[gid]]), CARHR2NAME[[gid]], "")
        gdesc <- ifelse(!is.null(CARHR2DESC[[gid]]), CARHR2DESC[[gid]], "")
    }
    theta.AR <- A / (A + R)
    theta.ar <- a / (a + r)
    log10fpkm.AR <- log10(A * (1/3) + R * (2/3) + 1)
    log10fpkm.ar <- log10(a + r + 1)
    df <- rbind(data.frame(A = log10fpkm.AR * cos(pi * theta.AR), R = log10fpkm.AR * sin(pi * theta.AR),
                           time = gsub(" hr", "", TIMELABELS), allopolyploid = "synthetic"),
                data.frame(A = log10fpkm.ar * cos(pi * theta.ar), R = log10fpkm.ar * sin(pi * theta.ar),
                           time = gsub(" hr", "", TIMELABELS), allopolyploid = "C. insueta"))
    g <- ggplot(df, aes(x = A, y = R, label = time, colour = allopolyploid, group = allopolyploid))
    g <- g + geom_abline(intercept = 0, slope = sin(pi/3)/cos(pi/3), colour = "#AAAAAA", linetype = 2)
    g <- g + geom_path() + theme_bw() + coord_fixed()
    g <- g + ggtitle(bquote(atop(.(paste0(gid, " ", gname)), atop(italic(.(gdesc)), "")))) 
    g <- g + geom_text(size = 4, colour = c(rep("#E41A1C", 9), rep("#377E88", 9))) 
    g <- g + xlab("log10(FPKM + 1) * cos(A-origin ratio * pi)") 
    g <- g + ylab("log10(FPKM + 1) * sin(A-origin ratio * pi)")
    g <- g + scale_color_manual(values = c("#B3CDE3", "#FBB4AE"))
    g <- g + xlim(min(c(-3, df$A, df$R), na.rm = T), max(c(3, df$A, df$R), na.rm = T)) + ylim(0, max(c(3, df$A, df$R), na.rm = T))
    g
}
    
PlotFPKMScatter.plot_rhxy <- function(A = NULL, R = NULL, a = NULL, r = NULL, gid, gname, gdesc) {
    if (is.null(A) && is.null(R)) {
        A = ca$fpkm$AA[gid, ] * (1/3)
        R = ca$fpkm$RR[gid, ] * (2/3)
        a = ca$fpkm$IA[gid, ]
        r = ca$fpkm$IR[gid, ]
        gname <- ifelse(!is.null(CARHR2NAME[[gid]]), CARHR2NAME[[gid]], "")
        gdesc <- ifelse(!is.null(CARHR2DESC[[gid]]), CARHR2DESC[[gid]], "")
    }
    theta.AR <- A / (A + R)
    theta.ar <- a / (a + r)
    log10fpkm.AR <- log10(A * (1/3) + R * (2/3) + 1)
    log10fpkm.ar <- log10(a + r + 1)
    #plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, max(c(3, log10fpkm.AR, log10fpkm.ar))),
    #     main = paste0(gid, " ", gname), sub = gdesc, xlab = "A-origin ratio", ylab = "log10(FPKM+1)")
    #abline(v = 1/3, col = "#AAAAAA", lty = 2)
    #lines(theta.AR, log10fpkm.AR, col = "red")
    #lines(theta.ar, log10fpkm.ar, col = "blue")
    #text(theta.AR, log10fpkm.AR, labels = gsub(" hr", "", TIMELABELS), col = "red")
    #text(theta.ar, log10fpkm.ar, labels = gsub(" hr", "", TIMELABELS), col = "blue")
    #legend("topleft", col = c("red", "blue"), legend = c("synthetic", "C. insueta"), lty = 1)
    df <- rbind(data.frame(AOriginRatio = theta.AR, log10FPKM = log10fpkm.AR, 
                           time = gsub(" hr", "", TIMELABELS), allopolyploid = "synthetic"),
                data.frame(AOriginRatio = theta.ar, log10FPKM = log10fpkm.ar,
                           time = gsub(" hr", "", TIMELABELS), allopolyploid = "C. insueta"))
    g <- ggplot(df, aes(x = AOriginRatio, y = log10FPKM, label = time, colour = allopolyploid, group = allopolyploid))
    g <- g + geom_vline(xintercept = 1/3, colour = "#AAAAAA", linetype = 2)
    g <- g + geom_path() + theme_bw()
    g <- g + ggtitle(bquote(atop(.(paste0(gid, " ", gname)), atop(italic(.(gdesc)), "")))) 
    g <- g + geom_text(size = 4, colour = c(rep("#E41A1C", 9), rep("#377E88", 9))) 
    g <- g + xlab("A-origin ratio") 
    g <- g + ylab("log10(FPKM + 1)")
    g <- g + scale_color_manual(values = c("#B3CDE3", "#FBB4AE"))
    g <- g + xlim(0, 1) + ylim(0, max(c(3, log10fpkm.AR, log10fpkm.ar), na.rm = T))
    g
}
    
    
PlotFPKMScatter <- function(ca) {
    deg.AA <- union(unlist(ca$DEG$UP$AA), unlist(ca$DEG$DW$AA))
    deg.RR <- union(unlist(ca$DEG$UP$RR), unlist(ca$DEG$DW$RR))
    deg.IA <- union(unlist(ca$DEG$UP$IA), unlist(ca$DEG$DW$IA))
    deg.IR <- union(unlist(ca$DEG$UP$IR), unlist(ca$DEG$DW$IR))
    deg.AL <- unique(c(deg.AA, deg.RR, deg.IA, deg.IR))
    for (g in rownames(ca$fpkm$AA)) {
        g.name <- ifelse(!is.null(CARHR2NAME[[g]]), CARHR2NAME[[g]], "")
        g.desc <- ifelse(!is.null(CARHR2DESC[[g]]), CARHR2DESC[[g]], "")
        if (length(grep(g, deg.AL)) > 0) {
            imgdir.xy <- paste0("~/Desktop/fpkm_scatter/xy/DEGs/", g, ".png")
            imgdir.pl <- paste0("~/Desktop/fpkm_scatter/pl/DEGs/", g, ".png")
            imgdir.rhxy <- paste0("~/Desktop/fpkm_scatter/rh/DEGs/", g, ".png")
        } else {
            imgdir.xy <- paste0("~/Desktop/fpkm_scatter/xy/nonDEGs/", g, ".png")
            imgdir.pl <- paste0("~/Desktop/fpkm_scatter/pl/nonDEGs/", g, ".png")
            imgdir.rhxy <- paste0("~/Desktop/fpkm_scatter/rh/nonDEGs/", g, ".png")
        }
        png(imgdir.rhxy, width = 400, height = 400)
        .p <- PlotFPKMScatter.plot_rhxy(A = ca$fpkm$AA[g, ] * (1/3), R = ca$fpkm$RR[g, ] * (2/3),
                         a = ca$fpkm$IA[g, ], r = ca$fpkm$IR[g, ], g, g.name, g.desc)
        plot(.p)
        dev.off()
        #png(imgdir.xy, width = 400, height = 400)
        #.p <- PlotFPKMScatter.plot_xy(A = ca$fpkm$AA[g, ] * (1/3), R = ca$fpkm$RR[g, ] * (2/3),
        #               a = ca$fpkm$IA[g, ], r = ca$fpkm$IR[g, ], g, g.name, g.desc)
        #plot(.p)
        #dev.off()
        #png(imgdir.pl, width = 530, height = 280)
        #.p <- PlotFPKMScatter.plot_rh(A = ca$fpkm$AA[g, ] * (1/3), R = ca$fpkm$RR[g, ] * (2/3),
        #               a = ca$fpkm$IA[g, ], r = ca$fpkm$IR[g, ], g, g.name, g.desc)
        #plot(.p)
        #dev.off()
    }
    
    
}


PlotRsqured <- function(ca) {
    rsq.AA_IA <- rsq.RR_IR <- rsq.AA_RR <- rsq.IA_IR <- rep(0, 9)
    
    log10A <- log10(ca$fpkm$AA[ca$EXP$genes$ALL, ] + 1)
    log10R <- log10(ca$fpkm$RR[ca$EXP$genes$ALL, ] + 1)
    log10a <- log10(ca$fpkm$IA[ca$EXP$genes$ALL, ] + 1)
    log10r <- log10(ca$fpkm$IR[ca$EXP$genes$ALL, ] + 1)
    
    for (i in 1:9) {
        rsq.AA_IA[i] <- summary(lm(log10A[, i] ~ log10a[, i]))$r.squared
        rsq.RR_IR[i] <- summary(lm(log10R[, i] ~ log10r[, i]))$r.squared
        rsq.AA_RR[i] <- summary(lm(log10A[, i] ~ log10R[, i]))$r.squared
        rsq.IA_IR[i] <- summary(lm(log10a[, i] ~ log10r[, i]))$r.squared
    }
    
    rsq <- rbind(rsq.AA_IA, rsq.RR_IR, rsq.AA_RR, rsq.IA_IR)
    rownames(rsq) <- c("AA - IA", "RR - IR", "AA - RR", "IA - IR")
    colnames(rsq) <- TIMELABELS
    rsq.df <- melt(rsq)
    colnames(rsq.df) <- c("combination", "time", "RSquared")
    gg <- ggplot(rsq.df, aes(x = time, colour = combination, y = RSquared, group = combination))
    gg <- gg + geom_line()
    gg <- gg + xlab("Time point") + ylab("R-squared value")
    gg <- gg + theme_bw()
    gg <- gg + scale_colour_brewer(palette = "Dark2")
    ca$EXP$rsquared.plot <- gg
    ca
}


GOHyper <- function(ca) {
    go <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")
    goid <- go[, 1]
    goid2genes <- apply(go, 1, function(x) { unlist(strsplit(x[2], ";")) })
    names(goid2genes) <- goid
    
    all.list <- ca$gene.names
    
    
    
    go <- list(UP = ca$DEG$UP, DW = ca$DEG$DW)
    
    for (sp in names(ca$DEG)) {
        for (ud in c("UP", "DW")) {
            for (tm in names(ca$DEG[[sp]][[ud]])) {
                deg.list <- ca$DEG[[sp]][[ud]][[tm]]
                for (goterm in names(goid2genes)) {
                    genes.in.goterm <- goid2genes[[goterm]]
                    has.goterm.null <- intersect(all.list, genes.in.goterm)
                    unhas.goterm.null <- setdiff(all.list, genes.in.goterm)
                    
                }
                
                
            }
        }
        
    }
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
    g <- g + facet_wrap( ~ time, ncol = 3)
    g <- g + xlab("A-origin ratio") + ylab("frequency")
    ca$BIAS$aorigin.ratio.dist <- g
    ca
}



EstBiasedHomeologs <- function(ca) {
    source(paste0(PATH$lib, '/homeoBiasGmm.R'))
    target.homeologs <- (ca$fpkm$IA + ca$fpkm$IR > FPKM_CUTOFF)
    bmat <- matrix("", ncol = 9, nrow = nrow(target.homeologs))
    colnames(bmat) <- TIMELABELS
    rownames(bmat) <- rownames(target.homeologs)
    
    for (i in 1:9) {
        d <- ca$aorigin.ratio[, i]
        d <- d[target.homeologs[, i]]
        d <- as.matrix(d[!is.na(d)])
        d[d == 0, 1] <- 1e-5
        d[d == 1, 1] <- 1 - 1e-5
        y <- d[, 1]
        
        is.group1 <- (y <= 0.05 | y > 0.95)
        y1 <- y[is.group1]
        y2 <- y[!is.group1]
        r <- c(length(y1)/length(y), length(y2)/length(y))
        p1 <- getParamMV(y1)
        p2 <- getParamMV(y2)
        p <- matrix(c(p1[1], p2[1], p1[2], p2[2]), 2)
        res <- fitBMM(y, r, p, tol=1E-15)
        cat('\n')
        
        pdfname <- paste0(path.expbias, '/hist_hratio_bmm_Cinsueta_', i, 'th_timepoints.pdf')
        plotBMM(y, res$ratio, res$param.beta, hist.by=0.02, my.xlab='A-origin ratio', pdf.name=pdfname)
        dev.off()
        cat('H-ratio (LR=1)\n')
        cat(paste0(getValsGivenLR(res$ratio, res$param.beta, lr=1), "\n"))
        loglike <- getLogLR(d, res$ratio, res$param.beta)
        bias.genes <- getBiasGenes(d, loglike$logLR, lr=1)
        
        bmat[target.homeologs[, i], i] <- 'T'
        bmat[rownames(bias.genes$df.bias.low), i] <- 'R'
        bmat[rownames(bias.genes$df.bias.high), i] <- 'A'
    }
    
    df.sk <- NULL
    u <- nrow(bmat[rowSums(bmat == 'A' | bmat == 'R') == 0, ] )
    for (i in 1:8) {
        is.from.A <- bmat[, i] == 'A'
        is.from.R <- bmat[, i] == 'R'
        is.from.N <- bmat[, i] != 'A' & bmat[, i] != 'R'
        is.to.A <- bmat[, i + 1] == 'A'
        is.to.R <- bmat[, i + 1] == 'R'
        is.to.N <- bmat[, i + 1] != 'A' & bmat[, i + 1] != 'R'
        tf <- TIMELABELS[i]
        tt <- TIMELABELS[i + 1]
        .df.sk <- data.frame(from = rep(c(paste0("A-biased ", tf), paste0("R-biased ", tf), paste0("non-biased ", tf)), each = 3)[-9],
                             to   = rep(c(paste0("A-biased ", tt), paste0("R-biased ", tt), paste0("non-biased ", tt)), times = 3)[-9],
                             weight = c(sum(is.from.A & is.to.A), sum(is.from.A & is.to.R), sum(is.from.A & is.to.N),
                                        sum(is.from.R & is.to.A), sum(is.from.R & is.to.R), sum(is.from.R & is.to.N),
                                        sum(is.from.N & is.to.A), sum(is.from.N & is.to.R)))
        df.sk <- rbind(df.sk, .df.sk)
    }
    skp <- gvisSankey(df.sk, from = "from", to = "to", weight = "weight",
                     options = list(width = "automatic", height = "automatic"))
    
    dfnum <- melt(apply(bmat, 2, table)[2:3, ])
    dfnum$Var1 <- ifelse(dfnum$Var1 == 'A', "A-biased", "R-biased")
    colnames(dfnum) <- c("homeolog", "time", "value")
    g <- ggplot(dfnum, aes(x = time, y = value, fill = homeolog, group = homeolog))
    g <- g + geom_bar(position = "dodge", stat = "identity")
    g <- g + ylab("Number of biased homeologs") + xlab("time")
    g <- g + theme_bw()
    
    
    ca$BIAS$biasedgenes.skplot <- skp
    ca$BIAS$num.of.homeologs <- g
    ca$BIAS$mat <- bmat
    
    write.table(apply(bmat, 2, table), file = paste0(path.expbias, "/number_of_biased_homeologs.xls"),
                row.names = T, col.names = T, sep = "\t", quote = FALSE)
    
    ca
}


BIAS.GO <- function(ca, p.cutoff = 1, q.cutoff = 1) {
    go <- list(ABiased = vector("list", 9), RBiased = vector("list", 9))
    names(go$ABiased) <- names(go$RBiased) <- TIMELABELS
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    for (i in 1:9) {
        .biased.a <- rownames(ca$BIAS$mat)[ca$BIAS$mat[, i] == "A"]
        .biased.r <- rownames(ca$BIAS$mat)[ca$BIAS$mat[, i] == "R"]
        .go.a <- DoGO(sig.genes = .biased.a, all.genes = ca$EXP$genes$ALL, p.cutoff = p.cutoff, q.cutoff = q.cutoff, ontology = "BP")
        .go.r <- DoGO(sig.genes = .biased.r, all.genes = ca$EXP$genes$ALL, p.cutoff = p.cutoff, q.cutoff = q.cutoff, ontology = "BP")
        .go.a <- .go.a$TABLE[ intersect(.go.a$TABLE$ID, target.terms), ]
        .go.r <- .go.r$TABLE[ intersect(.go.r$TABLE$ID, target.terms), ]
        .go.a$qvalue <- p.adjust(.go.a$pvalue, method = "BH")
        .go.r$qvalue <- p.adjust(.go.r$pvalue, method = "BH")
        go$ABiased[[i]] <- .go.a[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "qvalue", "geneID", "Count")]
        go$RBiased[[i]] <- .go.r[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "qvalue", "geneID", "Count")]
    }
    
    SaveExcel(go$ABiased, file.name = paste0(path.expbias, "/BIAS_GO_A-BiasedHomeologs.xlsx"))
    gc()
    SaveExcel(go$RBiased, file.name = paste0(path.expbias, "/BIAS_GO_R-BiasedHomeologs.xlsx"))
    gc()
    
    ca <- ca$BIASGO <- go
    ca
}

TAIRSTR2CARHR <- function(x) {
    unlist(TAIR2CARHR[unlist(strsplit(x, '/'))])
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
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "row", Colv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNCC2(), key.xlab = "log2(FPKM + 1)", margins = m,
              colsep = c(9, 18, 27, 36), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F,
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNCC2(), key.xlab = "log2(FPKM + 1)", margins = m,
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
    colnames(g.fpkm) <- paste0(rep(SPLABELS, each = 9), " (", TIMELABELS, ")")
    rownames(g.fpkm) <- df[, 3]
    g.fpkm <- g.fpkm[, -c(1, 9+1, 18+1, 27+1)]
    if (cls) {
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "row", Colv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    } else {
        heatmap.2(g.fpkm, trace = "none", scale = "none", dendrogram = "none", Colv = F, Rowv = F, key.title = "",
              hclustfun = function(x) hclust(x,method = 'ward.D2'),
              density.info = "none", col = COLSFUNC2(), key.xlab = "log2-foldchange", margins = m, keysize = 1.2,
              colsep = c(8, 16, 24, 32), sepwidth = c(0.1, 0.1), sepcolor = "white", cexRow = 0.8, cexCol = 0.8)
    }
}


CalcExpressedGenes <- function(ca) {
    exp.AA <- apply(ca$fpkm$AA, 2, function(x) {names(x)[x > FPKM_CUTOFF]})
    exp.RR <- apply(ca$fpkm$RR, 2, function(x) {names(x)[x > FPKM_CUTOFF]})
    exp.II <- apply(ca$fpkm$IA + ca$fpkm$IR, 2, function(x) {names(x)[x > FPKM_CUTOFF]})
    exp.AA.union <- unique(unlist(exp.AA))
    exp.RR.union <- unique(unlist(exp.RR))
    exp.II.union <- unique(unlist(exp.II))
    ca$EXP$genes <- list(AA = exp.AA.union, RR = exp.RR.union, II = exp.II.union,
                         ALL = unique(c(exp.AA.union, exp.AA.union, exp.II.union)))
    ca
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
                    ont = ontology, pAdjustMethod = "BH", keytype = "TAIR",
                    pvalueCutoff = p.cutoff, qvalueCutoff = q.cutoff, readable = FALSE)
        rgo <- summary(ego)
    list(TABLE = rgo, OBJ = ego)
}



DEG.GO <- function(ca, do.simplify = TRUE, cf = 0.7, p.cutoff = 1, q.cutoff = 1) {
    goobj <- go <- list(UP = ca$DEG$UP, DW = ca$DEG$DW)
    for (ud in c("UP", "DW")) {
        for (sp in names(ca$DEG[[ud]])) {
            for (tm in names(ca$DEG[[ud]][[sp]])) {
                .go <- DoGO(sig.genes = ca$DEG[[ud]][[sp]][[tm]], all.genes = ca$EXP$genes$ALL,
                            p.cutoff = p.cutoff, q.cutoff = q.cutoff, ontology = "BP")
                go[[ud]][[sp]][[tm]]    <- .go$TABLE
                goobj[[ud]][[sp]][[tm]] <- .go$OBJ
            }
        }
    }
    ca$DEGGO <- list(table = go, obj = goobj)
    
    target.terms <- read.table(paste0(PATH$dat, "/go/go2carhr.tsv"), header = F, sep = "\t")[, 1]
    gosim <- go
    for (ud in c("UP", "DW")) {
        for (sp in names(go[[ud]])) {
            for (tm in names(go[[ud]][[sp]])) {
                .tb <- go[[ud]][[sp]][[tm]]
                .target.terms <- intersect(.tb$ID, target.terms)
                .tb <- .tb[.target.terms,]
                .tb$qvalue <- p.adjust(.tb$pvalue, method = "BH")
                gosim[[ud]][[sp]][[tm]] <- .tb[, c("ID", "Description", "GeneRatio", "BgRatio",
                                                   "pvalue", "qvalue", "geneID", "Count")]
            }
        }
    }
    ca$DEGGOSMPL <- gosim
    
    for (ud in c("UP", "DW")) {
        for (sp in names(go[[ud]])) {
            SaveExcel(gosim[[ud]][[sp]], file.name = paste0(path.de, "/DEGGOSMPL_", ud, "_", sp, ".xlsx"))
            gc()
        }
    }
    
    
    
    
    ca
}




BindTairName <- function(dat) {
    if (!is.matrix(dat) && !is.data.frame(dat)) {
        dat <- as.matrix(dat)
        colnames(dat) <- "CARHR"
        rownames(dat) <- dat
    }
    rnm <- colnames(dat)
    tairid <- as.character(unlist(CARHR2TAIR[as.character(rownames(dat))]))
    genename <- as.character(unlist(CARHR2NAME[as.character(rownames(dat))]))
    genedesc <- as.character(unlist(CARHR2DESC[as.character(rownames(dat))]))
    dat <- data.frame(dat, TAIR = tairid, name = genename, description = genedesc, stringsAsFactors = F)
    dat[is.na(dat)] <- ""
    colnames(dat) <- c(rnm, "TAIR", "name", "description")
    dat
}








TestRatioShifts <- function(ca) {
    pval <- matrix(NA, ncol = 6, nrow = nrow(ca$fpkm$AA))
    rownames(pval) <- rownames(ca$fpkm$AA)
    colnames(pval) <- c("A_is_0", "R_is_0", "a_is_0", "r_is_0", "A_is_a", "R_is_r")
    A <- log10(ca$fpkm$AA * (1/3) + 1)
    R <- log10(ca$fpkm$RR * (2/3) + 1)
    a <- log10(ca$fpkm$IA + 1)
    r <- log10(ca$fpkm$IR + 1)
    
    for (g in rownames(ca$fpkm$AA)) {
        try.A_is_0 <- try(t.test(A[g, ])$p.value)
        try.R_is_0 <- try(t.test(R[g, ])$p.value)
        try.a_is_0 <- try(t.test(a[g, ])$p.value)
        try.r_is_0 <- try(t.test(r[g, ])$p.value)
        try.A_is_a <- try(t.test(A[g, ], a[g, ])$p.value)
        try.R_is_r <- try(t.test(R[g, ], r[g, ])$p.value)
        if (class(try.A_is_0) != "try-error") {pval[g, 1] <- try.A_is_0}
        if (class(try.R_is_0) != "try-error") {pval[g, 2] <- try.R_is_0}
        if (class(try.a_is_0) != "try-error") {pval[g, 3] <- try.a_is_0}
        if (class(try.r_is_0) != "try-error") {pval[g, 4] <- try.r_is_0}
        if (class(try.A_is_a) != "try-error") {pval[g, 5] <- try.A_is_a}
        if (class(try.R_is_r) != "try-error") {pval[g, 6] <- try.R_is_r}
    }
    
    qval <- apply(pval, 2, p.adjust, method = "BH")
    pval[is.na(pval)] <- 1
    qval[is.na(qval)] <- 1
    
    .cutoff.q <- 0.001 
    A_isnot_0 <- (qval[, "A_is_0"] < .cutoff.q)
    R_isnot_0 <- (qval[, "R_is_0"] < .cutoff.q)
    a_isnot_0 <- (qval[, "a_is_0"] < .cutoff.q)
    r_isnot_0 <- (qval[, "r_is_0"] < .cutoff.q)
    A_isnot_a <- (qval[, "A_is_a"] < .cutoff.q)
    R_isnot_r <- (qval[, "R_is_r"] < .cutoff.q)
    
    
    left_to_right <- qval[(!R_isnot_r) & (A_isnot_a) & (!a_isnot_0) & (A_isnot_0), ]
    right_to_left <- qval[(!R_isnot_r) & (A_isnot_a) & (a_isnot_0) & (!A_isnot_0), ]
    bottom_to_top <- qval[(!A_isnot_a) & (R_isnot_r) & (!r_isnot_0) & (R_isnot_0), ]
    top_to_bottom <- qval[(!A_isnot_a) & (R_isnot_r) & (r_isnot_0) & (!R_isnot_0), ]
    
    for (gid in rownames(left_to_right)) {
        png(paste0("~/Desktop/fpkm_scatter/left_to_right/", gid, ".xy.png"), 400, 400)
        plot(PlotFPKMScatter.plot_xy(gid = gid))
        dev.off()
        png(paste0("~/Desktop/fpkm_scatter/left_to_right/", gid, ".rh.png"), 400, 400)
        plot(PlotFPKMScatter.plot_rhxy(gid = gid))
        dev.off()
    }
    for (gid in rownames(right_to_left)) {
        png(paste0("~/Desktop/fpkm_scatter/right_to_left/", gid, ".xy.png"), 400, 400)
        plot(PlotFPKMScatter.plot_xy(gid = gid))
        dev.off()
        png(paste0("~/Desktop/fpkm_scatter/right_to_left/", gid, ".rh.png"), 400, 400)
        plot(PlotFPKMScatter.plot_rhxy(gid = gid))
        dev.off()
    }
    for (gid in rownames(bottom_to_top)) {
        png(paste0("~/Desktop/fpkm_scatter/bottom_to_top/", gid, ".xy.png"), 400, 400)
        plot(PlotFPKMScatter.plot_xy(gid = gid))
        dev.off()
        png(paste0("~/Desktop/fpkm_scatter/bottom_to_top/", gid, ".rh.png"), 400, 400)
        plot(PlotFPKMScatter.plot_rhxy(gid = gid))
        dev.off()
    }
    for (gid in rownames(top_to_bottom)) {
        png(paste0("~/Desktop/fpkm_scatter/top_to_bottom/", gid, ".xy.png"), 400, 400)
        plot(PlotFPKMScatter.plot_xy(gid = gid))
        dev.off()
        png(paste0("~/Desktop/fpkm_scatter/top_to_bottom/", gid, ".rh.png"), 400, 400)
        plot(PlotFPKMScatter.plot_rhxy(gid = gid))
        dev.off()
    }
    
    
    
    
}





PlotDEGVenn <- function(ca) {
    for (ud in c("UP", "DW")) {
        for (tm in 1:8) {
            .dat.A <- ca$DEG[[ud]]$AA[[tm]]
            .dat.a <- ca$DEG[[ud]]$IA[[tm]]
            .dat.r <- ca$DEG[[ud]]$IR[[tm]]
            .dat.R <- ca$DEG[[ud]]$RR[[tm]]
            .dat <- rbind(data.frame(elements = .dat.A, sets  = 'AA'),
                          data.frame(elements = .dat.a, sets  = 'Ia'),
                          data.frame(elements = .dat.r, sets  = 'Ir'),
                          data.frame(elements = .dat.R, sets  = 'RR'))
            png(paste0(path.de, "/DEG_VENN1_", ud, "_", tm, "th-timepoints.png"))
            plot(venn(list(AA = .dat.A, IA = .dat.a, IR = .dat.r, RR = .dat.R)))
            dev.off()
            png(paste0(path.de, "/DEG_VENN2_", ud, "_", tm, "th-timepoints.png"))
            plot(venneuler(.dat))
            dev.off()
        }
    }
}



PlotAll <- function(ca) {
    ca <- PlotAoriginRatioDist(ca)
    ca <- IdentifyDEGsFCPlot(ca, q.cutoff = 0.1)
    ca <- PlotRsqured(ca)
    
    
}












.common.RData  <- paste0(PATH$lib, "/src/common.RData")
.ca.RData <- paste0(PATH$lib, "/src/ca.RData")



stop()
if (file.exists(.ca.RData)) {
    load(.common.RData)
    load(.ca.RData)
} else {
    design <- GetDesign()
    ca <- NULL
    ca <- GetCounts(ca, design)
    ca <- CalcFPKM(ca)
    ca <- CalcExpressedGenes(ca)
    ca <- EstBiasedHomeologs(ca)
    ca <- IdentifyDEGs(ca, q.cutoff = 0.1)
    ca <- DEG.GO(ca)
    ca <- BIAS.GO(ca)
    
    ca <- PlotAll(ca)
 
    # PlotFPKMScatter(ca)  # run this if needed
    
    
   
    save(ca, design, file = .ca.RData)
}





















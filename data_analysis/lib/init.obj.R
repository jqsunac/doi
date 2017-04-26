## 
## Create hash table for ID mapping.
## Only run this file for fisrt time, otherwise please load the .RData.
## It takes more about 15 minutes to run this file.
##
message("Start to run common.initobj.R.")

.f.id <- paste0(PATH$lib, "/src/chi_m25.txt")


## ./get_GO.py > cardamine.go.txt
## awk 'BEGIN{FS="\t"}{printf "\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\n", $1, $2, $5, $3, $4, $7}' chi_m25.txt > cardamine.id.txt
if (! file.exists(.f.id)) stop("Error: Use akw to create basis file!")
#if (! file.exists(.f.go)) stop("Error: Run get_GO.py!")

library(KEGG.db)
library(org.At.tair.db)
ls("package:org.At.tair.db")


#--------------------------------------------------------------
# Gene ID hash table
#--------------------------------------------------------------

# CARHR2TAIR; CARHR2NAME; CARHR2DESC;
#--------------------------------------------------------------
message("Create CARHR2TAIR, CARHR2NAME, and CARHR2DESC.")
carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
carhrid <- unique(as.character(carhrf[, 1]))
CARHR2TAIR <- vector("list", length(carhrid))
CARHR2NAME <- vector("list", length(carhrid))
CARHR2DESC <- vector("list", length(carhrid))
names(CARHR2TAIR) <- names(CARHR2NAME) <- names(CARHR2DESC) <- carhrid
for (i in 1:length(CARHR2TAIR)) {
  CARHR2TAIR[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 2])
  CARHR2NAME[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 5])
  CARHR2DESC[[i]] <- as.character(carhrf[carhrf[, 1] == carhrid[i], 6])
}

# TAIR2CARHR
#--------------------------------------------------------------
message("Creating TAIR2CARHR ...")
carhrf <- read.table(.f.id, header = TRUE, sep = "\t", fill = TRUE)
carhrf[, 6] <- gsub(" \\[.+\\]$", "", carhrf[,6])
tairid <- unique(as.character(carhrf[, 2]))
TAIR2CARHR <- vector("list", length(tairid))
names(TAIR2CARHR) <- tairid
for (i in 1:length(TAIR2CARHR)) {
  TAIR2CARHR[[i]] <- as.character(carhrf[carhrf[, 2] == tairid[i], 1])
}


CREATED_DATE <- date()
SESS_INFO <- sessionInfo()

save(CARHR2TAIR, TAIR2CARHR, CARHR2NAME, CARHR2DESC, CREATED_DATE, SESS_INFO,
     file = .common.RData)
# till here



#--------------------------------------------------------------
# KEGG hash table
#--------------------------------------------------------------

# KEGG2TAIR and TAIR2KEGG
#--------------------------------------------------------------
message("Creating TAIR2KEGG ...")
KEGG2TAIR <- as.list(org.At.tairARACYC)

message("Creating KEGG2TAIR ...")
.kegg_names <- names(KEGG2TAIR)
.tair_names <- unique(as.character(unlist(KEGG2TAIR)))
TAIR2KEGG <- vector("list", length(.tair_names))
names(TAIR2KEGG) <- .tair_names
for (i in 1:length(KEGG2TAIR)) {
    TAIR2KEGG[KEGG2TAIR[[i]]] <- lapply(TAIR2KEGG[KEGG2TAIR[[i]]], function(x) {c(x, .kegg_names[i])})
}

# CARHR2KEGG and KEGG2CARHR
#--------------------------------------------------------------
message("Creating CARHR2KEGG ...")
tairid <- unique(names(TAIR2KEGG))
carhrid <- as.character(unlist(TAIR2CARHR[tairid]))
CARHR2KEGG <- vector("list", length(carhrid))
names(CARHR2KEGG) <- carhrid
for (i in 1:length(CARHR2KEGG)) {
  CARHR2KEGG[[i]] <- as.character(TAIR2KEGG[[as.character(CARHR2TAIR[carhrid[i]])]])
}

message("Creating KEGG2CARHR ...")
kgid <- names(KEGG2TAIR)
KEGG2CARHR <- vector("list", length(kgid))
names(KEGG2CARHR) <- kgid
for (i in 1:length(kgid)) {
  KEGG2CARHR[[i]] <- as.character(unlist(TAIR2CARHR[KEGG2TAIR[[kgid[i]]]]))
}



#--------------------------------------------------------------
# GO hash table
#--------------------------------------------------------------

# TAIR2GO and GO2TAIR
#--------------------------------------------------------------
message("Creating TAIR2GO ...")
.go <- as.list(org.At.tairGO)
.tair_names <- names(.go)
.go_names   <- unique(as.character(unlist(lapply(.go, names))))
TAIR2GO <- vector("list", length(.tair_names))
names(TAIR2GO) <- .tair_names
for (i in 1:length(TAIR2GO)) {
  TAIR2GO[[i]] <- names(.go[[names(TAIR2GO)[i]]])
}

message("Creating GO2TAIR ...")
GO2TAIR <- as.list(org.At.tairGO2ALLTAIRS)


#--------------------------------------------------------------
# GO2CARHR and CARHR2GO
#--------------------------------------------------------------
message("Creating CARHR2GO ...")
tairid <- unique(names(TAIR2GO))
carhrid <- as.character(unlist(TAIR2CARHR[tairid]))
CARHR2GO <- vector("list", length(carhrid))
names(CARHR2GO) <- carhrid
for (i in 1:length(CARHR2GO)) {
  CARHR2GO[[i]] <- as.character(TAIR2GO[[as.character(CARHR2TAIR[carhrid[i]])]])
}

message("Creating GO2CARHR ...")
kgid <- names(GO2TAIR)
GO2CARHR <- vector("list", length(kgid))
names(GO2CARHR) <- kgid
for (i in 1:length(kgid)) {
  GO2CARHR[[i]] <- as.character(unlist(TAIR2CARHR[GO2TAIR[[kgid[i]]]]))
}





#--------------------------------------------------------------
# KEGG
#--------------------------------------------------------------
KEGGNAME2ID <- as.list(KEGGPATHNAME2ID)
KEGGID2NAME <- as.list(KEGGPATHID2NAME)







message(paste0("Saving all hash table in to ", .common.RData))


TAIR2NAME, TAIR2DESC,
     TAIR2KEGG, KEGG2TAIR, CARHR2KEGG, KEGG2CARHR, KEGGNAME2ID, KEGGID2NAME,
     TAIR2GO, GO2TAIR, CARHR2GO, GO2CARHR, CREATED_DATE, SESS_INFO,
     file = .common.RData)



message("Finish running common.initobj.R.")


# The code below performs the two-dimensional QTL scan over the specified genomic region. 
# We analyze the chromosome 2 hotspot from Keller et al 2018 (GENETICS).
# Each scan considers two traits. 
##First read in the arguments listed at the command line
n_chr <- 7
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
run_num <- as.numeric(args$run_num)
print(run_num)
(nsnp <- as.numeric(args$nsnp))
(s1 <- as.numeric(args$s1))

###############


# load expression traits
readRDS(paste0("Chr", n_chr, "hot_local.rds")) -> local


nlocal <- ncol(local)
# create matrix of two expression traits
inds <- combn(nlocal, 2)[ , proc_num + 1]
pheno <- local[ , inds] #inds has length 2

# load chr2 allele probabilities
readRDS(paste0("Chr", n_chr, "_aprobs.rds")) -> geno # genoprobs_chr2.rds is on SQUID

# load kinship matrix (LOCO, ie, for chromosome 2, ie, doesn't use chr2 data)
readRDS(paste0("Chr", n_chr, "_kinship.rds")) -> kinship

# load covariates
readRDS("addcovar.rds") -> covar


# remove subjects with missing data

id2keep <- rownames(local)
gg <- geno
gg2 <- gg[rownames(gg) %in% id2keep, , ]
kk <- kinship
kk2 <- kk[rownames(kk) %in% id2keep, colnames(kk) %in% id2keep]
cc2 <- covar[rownames(covar) %in% id2keep, ]


# verify that names match in all objects
sum(rownames(pheno) == rownames(gg2))
sum(rownames(pheno) == rownames(kk2))
sum(rownames(pheno) == colnames(kk2))
sum(rownames(pheno) == rownames(cc2))
phenames <- c(colnames(local)[local_indic], colnames(hotspot)[hot_indic])
# two-dimensional scan

library(qtl2pleio)
s_out <- scan_pvl(probs = gg2,
         pheno = pheno,
         kinship = kk2,
         addcovar = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
         start_snp = s1,
         n_snp = nsnp
           )
# make a profile lod tibble
map <- readRDS("map.rds")
s_out <- tidy_scan_pvl(s_out, pmap = map[[as.character(n_chr)]]) # changes for each hotspot

# write output
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(s_out, fn_out, quote = FALSE)
q("no")

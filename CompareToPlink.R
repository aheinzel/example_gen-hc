library(dplyr)
library(magrittr)
library(purrr)
source("GenToHardCall.R")

system("plink --gen sample_data.gen --sample sample_data.samples --oxford-single-chr 1 --recodeA --out plink")
df.plink.raw <- read.table("plink.raw", sep = " ", header = T, check.names = F) %>%
  select(-c(FID, PAT, MAT, SEX, PHENOTYPE))

df.gen <- read.table("sample_data.gen", sep = " ")
df.samples <- read.table("sample_data.samples")
df.hc <- oxfordGenToHardCalls(df.gen, 0.1)

samples <- df.samples[-c(1:2), 2]
colnames(df.hc)[6:ncol(df.hc)] <- samples

rownames(df.hc) <- df.hc %$% paste(V2, V5, sep = "_")
df.hc <- df.hc[-c(1:5)]
df.hc <- t(df.hc)
df.hc <- cbind(data.frame(IID = rownames(df.hc)), df.hc)
rownames(df.hc) <- NULL


stopifnot(sub("_[^_]+$", "", colnames(df.hc)) == sub("_[^_]+$", "", colnames(df.plink.raw)))
df.hc.flip <- lapply(
  2:ncol(df.hc),
  function(x){
    if(colnames(df.hc)[x] == colnames(df.plink.raw)[x]){
      return(df.hc[[x]])
    }else{
      return(2 - df.hc[[x]])
    }
  }
)

names(df.hc.flip) <- sub("_[^_]+$", "_X", colnames(df.hc)[-1])
df.hc.flip <- cbind(data.frame(IID = df.hc$IID), do.call(cbind, df.hc.flip))

all((is.na(df.hc.flip) & is.na(df.plink.raw)) | df.hc.flip ==  df.plink.raw)

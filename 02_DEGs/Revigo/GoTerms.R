#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script created in version R 4.2.2 on xx/xx/2023
# This script:
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
setwd("C:/Users/NJCB/Documents/Purdue/daphnia_RNAseq/results/02_degs/Revigo_NO_CBSWP/")



list.files()

#dat <- read.csv("168.genotype_clone-temp.merged_annotations.tsv", sep="\t") 
dat <- read.csv("./168.NO_CBSWP.genotype_clone-temp.merged_annotations.tsv", sep="\t") 


#=============================================================================================================#


colnames(dat)

gos <- dat$GOs

OUT <- NULL
for(n in 1:length(gos)){
  go <- gos[n]
  go <- unlist(strsplit(go, ","))
  if(length(go) == 0) {next}
  out <- cbind(1:length(go), go)
  OUT <- rbind(OUT, out)
}

gogo <- OUT
gogo <- gogo[order(gogo[, 2]), ]
gogo <- gogo[-c(1:204), ]
head(gogo)

gogo <- data.frame(table(gogo[, 2]))
gogo <- gogo[order(gogo[, 2], decreasing = TRUE), ]

nrow(gogo[which(gogo[, 2] > 5), ]) # keep those with 5 or more counts
gogo2 <- gogo[which(gogo[, 2] > 5), ]
write.table(gogo2, "go_168_NO_CBSWP.txt", sep="\t", append = FALSE,col.names=FALSE, row.names = FALSE,  quote = FALSE)

#gogo2 <- gogo[which(gogo[, 2] > 1), ] # keep those with 2 or more counts
#write.table(gogo2, "go_168.txt", sep="\t", append = FALSE,col.names=FALSE, row.names = FALSE,  quote = FALSE)

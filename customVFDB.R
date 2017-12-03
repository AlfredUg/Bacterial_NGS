vfid <- read.table("vf_id.txt")
vf_acc <- read.table("vf_acc.txt", quote="\"", comment.char="")
vfshort <- read.table("vfshort.txt")
vf_full <- read.delim2("vf_full.txt", header=FALSE, quote="")
vf_pathogen <- read.delim("vf_pathogen.txt", header=FALSE)
vf_cat <- read.delim("vf_cat.txt", header=FALSE)

vfdb <- cbind(vfid, vf_acc, vf_full, vf_pathogen, vfshort, vf_cat)

names(vfdb) <- c("vfid","vfacc", "vf_full_name", "vf_pathogen", "vf_short_name", "vf_category")

write.csv(vfdb, file="vfdb.csv")

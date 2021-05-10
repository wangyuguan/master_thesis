C <- read.table("mat4chol.txt",header=F,sep=" ")
cmat <- chol(C)
cmat <- t(cmat)
write.table(cmat, file ="chol_res.txt", row.names=F, col.names=F)
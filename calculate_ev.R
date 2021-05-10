C <- read.table("corr_matrix.txt",header=F,sep=" ")
ev <- eigen(C)$values
write.table(ev, file ="ev_res.txt", row.names=F, col.names=F)

C <- read.table("mat4inv.txt",header=F,sep=" ")
mat_inv <- solve(C)
write.table(mat_inv, file ="inv_res.txt", row.names=F, col.names=F)
popInfo = read.csv("../external_data/data.clst", sep = "\t", header = FALSE, stringsAsFactors=FALSE)
#popInfo_new <- popInfo[!is.na(match(popInfo[,3],c("Yamnaya","LBK_EN","Loschbour"))),];
popFam = read.csv("../external_data/haak_fig3.LDprune.fam", sep = " ", header = FALSE, stringsAsFactors=FALSE )

gm_imputed <- data.frame(fread('../external_data/imputed_haak_fig3_genotypes.txt'),row.names=1);
out <- do.call(rbind,lapply(c("Yamnaya","LBK_EN","Loschbour"), function(x)
  {
    fids = popInfo$V2[which(popInfo[,3] == x)]
    inds = which(!is.na(match(popFam$V2,fids)))
    temp <- rbind(gm_imputed[inds,],gm_imputed[inds,]); ## replicate the rows to avoid having just 1 sample 
    freqMat = colMeans(temp)/2;
    return(freqMat)
}))

write.table(out, '../external_data/ancestral_freqs.txt');

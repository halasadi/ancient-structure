popInfo = read.csv("../external_data/data.clst", sep = "\t", header = FALSE, stringsAsFactors=FALSE)
#popInfo_new <- popInfo[!is.na(match(popInfo[,3],c("Yamnaya","LBK_EN","Loschbour"))),];
popFam = read.csv("../external_data/haak_fig3_set.fam", sep = " ", header = FALSE, stringsAsFactors=FALSE )

out <- do.call(rbind,lapply(c("Yamnaya","LBK_EN","Loschbour"), function(x)
  {
    fids = popInfo$V2[which(popInfo[,3] == x)]
    inds = which(!is.na(match(popFam$V2,fids)))
    temp <- rbind(gm[inds,],gm[inds,]); ## replicate the rows to avoid having just 1 sample 
    freqMat = colMeans(temp)/2;
    return(freqMat)
}))

write.table(out, '../external_data/ancestral_freqs.txt');

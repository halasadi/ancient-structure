library(snpStats)
reads <- read.plink(bed="../external_data/haak_fig3_set.bed",
                    bim="../external_data/haak_fig3_set.bim",
                    fam="../external_data/haak_fig3_set.fam")
gm = as(reads$genotypes, 'numeric')

gm_imputed <- apply(gm,2,function(x)
           {
              y=x;
              y[is.na(x)]=rbinom(1,2,mean(x,na.rm=TRUE)/2);
              return(y);
            })

write.table(gm_imputed, '../external_data/imputed_haak_fig3_genotypes.txt');



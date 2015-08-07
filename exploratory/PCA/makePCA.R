library(data.table)

# read in the data
data <- data.frame(fread("../external_data/imputed_haak_fig3_genotypes.txt"),row.names=TRUE);

pca = prcomp(data, scale = FALSE)


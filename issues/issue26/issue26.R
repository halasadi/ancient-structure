
## issue 26

setwd('~/Documents/ancient-structure/issues/issue26/')
library(data.table)
library(CountClust)
data.clst <- read.table('../../internal_data/indians_admixture/haak_data_indians_ee_ca.clst', stringsAsFactors = FALSE, header=FALSE, sep = " ");
data.fam <- read.table('../../internal_data/indians_admixture/haak_fig3_indians.LDprune.fam');

pop_ids <- data.clst[match(data.fam$V2, data.clst$V2),3]


plot_admixture <- function(K, pop_ids)
{
  P_data <- data.frame(fread(paste0('../../internal_data/indians_admixture/haak_fig3_indians.LDprune.',K,'.P')));
  Q_data <- data.frame(fread(paste0('../../internal_data/indians_admixture/haak_fig3_indians.LDprune.',K,'.Q')))
  
  modern_Q <- Q_data[1:620,];
  ancient_Q <- Q_data[-(1:620),];
  
  pop_ids_modern <- as.matrix(pop_ids[1:620]);
  
  
  if(!dir.exists('Admixture/modern')) dir.create('Admixture/modern')
  StructureObj_omega(modern_Q, samp_metadata = pop_ids_modern,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/modern',
                     control = list(cex.axis=1))
  
  if(!dir.exists('Admixture/ancient')) dir.create('Admixture/ancient')
  pop_ids_ancient <- as.matrix(pop_ids[-(1:620)]);
  StructureObj_omega(ancient_Q, samp_metadata = pop_ids_ancient,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/ancient',
                     control = list(cex.axis=1))

}

plot_admixture(3,pop_ids)
plot_admixture(5,pop_ids)
plot_admixture(7,pop_ids)

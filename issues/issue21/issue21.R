
## issue 21

setwd('~/Documents/ancient-structure/issues/issue21/')
library(data.table)
library(CountClust)
data.clst <- read.csv('../../external_data/data.clst', stringsAsFactors = FALSE, header=FALSE, sep = "\t");
data.fam <- read.table('../../external_data/haak_fig3.LDprune.fam');

pop_ids <- data.clst[match(data.fam$V2, data.clst$V2),3]

plot_admixture <- function(K, pop_ids)
{
  Q_data <- data.frame(fread('../../internal_data/YAM_LBK_GREEN.txt'))[,-1];
 
  modern_Q <- Q_data[1:331,];
  ancient_Q <- Q_data[-(1:331),];
  
  pop_ids_modern <- as.matrix(pop_ids[1:331]);
  
  
  if(!dir.exists('Admixture/modern_yam_lbk_green')) dir.create('Admixture/modern_yam_lbk_green')
  StructureObj_omega(modern_Q, samp_metadata = pop_ids_modern,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/modern_yam_lbk_green',
                     control = list(cex.axis=1))
  
  if(!dir.exists('Admixture/ancient_yam_lbk_green')) dir.create('Admixture/ancient_yam_lbk_green')
  pop_ids_ancient <- as.matrix(pop_ids[-(1:331)]);
  StructureObj_omega(ancient_Q, samp_metadata = pop_ids_ancient,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/ancient_yam_lbk_green',
                     control = list(cex.axis=1));
}

plot_admixture(3,pop_ids)

library(data.table)
library(CountClust)
data.clst <- read.csv('../external_data/data.clst', stringsAsFactors = FALSE, header=FALSE, sep = "\t");
data.fam <- read.table('../external_data/haak_fig3.LDprune.fam');

pop_ids <- data.clst[match(data.fam$V2, data.clst$V2),3]


plot_admixture <- function(K, pop_ids)
{
  #P_data <- data.frame(fread('Admixture/YAM_LBK_GREEN.txt'))
  Q_data <- read.table('Admixture/YAM_LBK_UNKNOWN.txt')
  
  modern_Q <- Q_data[1:331,];
  ancient_Q <- Q_data[-(1:331),];
  
  pop_ids_modern <- cbind.data.frame(pop_ids[1:331]);
  colnames(pop_ids_modern) <-"pop"
  
  if(!dir.exists('Admixture/modern_issue_21_unknown')) dir.create('Admixture/modern_issue_21_unknown')
  StructureObj_omega(modern_Q, samp_metadata = pop_ids_modern,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/modern_issue_21_unknown',
                     control = list(cex.axis=1))
  
  if(!dir.exists('Admixture/ancient_issue_21_unknown')) dir.create('Admixture/ancient_issue_21_unknown')
  pop_ids_ancient <- cbind.data.frame(pop_ids[-(1:331)]);
  colnames(pop_ids_ancient) <-"pop"
  StructureObj_omega(ancient_Q, samp_metadata = pop_ids_ancient,
                     partition = 'TRUE',batch_lab = NULL,
                     path_struct = 'Admixture/ancient_issue_21_unknown',
                     control = list(cex.axis=1))
}

plot_admixture(3,pop_ids)
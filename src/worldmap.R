## Ancient Structure on a world map

locs <- read.table("../external_data/locs.txt", sep=",", header=TRUE);
pop_ids <- as.vector(as.matrix(read.table("../external_data/pop_ids_Haak_et_al.txt")));
pop_ids_modern <- pop_ids[1:331];


Q_data <- read.table("../bin/admixture_modern/haak_fig3.LDprune.moderns.5.Q");

Q_modern <- Q_data;

Q.frame <- data.frame(pop_ids_modern, Q_modern);

library(dplyr)

mean_admix <- Q.frame %>% group_by(pop_ids_modern) %>% summarise_each(funs(mean))

mean_omega_mat <- as.matrix(mean_admix[,-1]);

match_indices <- match(locs[,1], mean_admix[,1]$pop_ids_modern)

mean_omega_mat <- mean_omega_mat[match_indices,];

omega.frame <- data.frame(locs,mean_omega_mat);

library(maps)
library(mapdata)
library(mapplots)
library(scales)

png(filename="../plots/geostructure_admixture_modern_5.png",res=200)
map("worldHires",
    ylim=c(35,70), xlim=c(-25,42), # Re-defines the latitude and longitude range
    col = "gray", fill=TRUE, mar=c(0.1,0.1,0.1,0.1))

lapply(1:dim(omega.frame)[1], function(r) 
  add.pie(z=as.integer(100*omega.frame[r,-(1:3)]), 
          x=omega.frame$long[r], y=omega.frame$lat[r], labels=c("","",""),
          col=c(alpha(2,0.6),alpha(3,0.6),alpha(4,0.6),alpha(5,0.6),alpha(6,0.6))));
dev.off()



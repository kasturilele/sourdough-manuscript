#script for all the figures + supplementary figures in the paper

#--------- figure 2: pairwise coexistence ------------

library(ggplot2)
library(dplyr)
library(cowplot)
#reading data from the pairwise parameters sheet - new data
paired_data <-  read.table(file = "/Users/kasturilele/Documents/community/est_all_pairs_6-3.txt", sep = ",", header = TRUE)

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum","C.paralimentarius","S.cerevisiae", "W.anomalus","K.humilis","K.servazzii")
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")
colour_strainID <- c("163" = "#88CCEE",
                     "228" = "#332288",
                     "253" = "#117733",
                     "177" = "#44AA99",
                     "17B2" = "#FDAE99",
                     "0092a" = "#882255",
                     "232" = "#DDCC77",
                     "460" = "#CC6677",
                     "550" = "#AA4499") #colour names with strain IDs

#predictions of coexistence
#making a few changes to separate bacteria from yeasts
bacteria <- c("17B2","0092a","232","460","550")
yeasts <- c("253","163","228","177")
paired_data <- subset(paired_data, Replicate < 20)

coexist <- data.frame(Strain_1=character(),
                      Strain_2=character(),
                      rho=double(),
                      f2_f1=double(),
                      b_y=character())

pred_coexist <- data.frame(Strain_1=character(),
                           Strain_2=character(),
                           coexist = integer(),
                           S1_win = integer(),
                           S2_win = integer(),
                           priority = integer())
c <- 1
for(i in 1:8){
  S1 <- strain_order[i]
  k <- i+1
  for(j in k:9){
    S2 <- strain_order[j]
    flip <- FALSE
    temp_pars <- subset(paired_data, Strain_1 == S1 & Strain_2 == S2)
    print(nrow(temp_pars))
    if(nrow(temp_pars) < 1){
      temp_pars <- subset(paired_data, Strain_1 == S2 & Strain_2 == S1)
      flip <- TRUE #S1 and S2 are flipped, this affects calculation of f2/f1
    }
    
    #create new blank data frame
    temp_coexist <- data.frame(Strain_1=character(),
                               Strain_2=character(),
                               rho=double(),
                               f2_f1=double())
    N <- nrow(temp_pars)
    for(r in 1:N){
      #calculate rho and f2/f1 for all pairs
      temp_rho <- (temp_pars[r,"a12"]/temp_pars[r,"a11"])*(temp_pars[r,"a21"]/temp_pars[r,"a22"])
      temp1_f2f1 <- (temp_pars[r,"a11"]/temp_pars[r,"a22"])*(temp_pars[r,"a12"]/temp_pars[r,"a21"])
      temp_f2f1 <- sqrt(abs(temp1_f2f1))*(temp_pars[r,"r2"]/temp_pars[r,"r1"])
      
      
      if(flip == TRUE){
        temp_coexist[r,1] <- temp_pars[r, "Strain_2"]
        temp_coexist[r,2] <- temp_pars[r, "Strain_1"]
        temp_coexist[r,3] <- sqrt(abs(temp_rho))
        temp_coexist[r,4] <- (1/temp_f2f1)
        #check for cross kingdom interactions
        if(temp_pars$Strain_2[1] %in% bacteria & temp_pars$Strain_1[1] %in% yeasts) {
          temp_coexist[r,5] <- "yes"
        }else{
          temp_coexist[r,5] <- "no"
        }
      }else{
        temp_coexist[r,1] <- temp_pars[r, "Strain_1"]
        temp_coexist[r,2] <- temp_pars[r, "Strain_2"]
        temp_coexist[r,3] <- sqrt(abs(temp_rho))
        temp_coexist[r,4] <- (temp_f2f1)
        #check for cross kingdom interactions
        if(temp_pars$Strain_1[1] %in% bacteria & temp_pars$Strain_2[1] %in% yeasts) {
          temp_coexist[r,5] <- "yes"
        }else{
          temp_coexist[r,5] <- "no"
        }
      }
    }
    
    #from the filled in data frame, calculate predictions of coexistence
    rho <- temp_coexist$rho
    f2_f1 <- temp_coexist$f2_f1
    var_coexist <- f2_f1 > rho & f2_f1 < 1/rho
    var_sp1 <- f2_f1 < rho & f2_f1 < 1/rho
    var_sp2 <- f2_f1 > rho & f2_f1 > 1/rho
    var_prior <- f2_f1 < rho & f2_f1 > 1/rho
    
    pred_coexist[c,1] <- S1
    pred_coexist[c,2] <- S2
    pred_coexist[c,3] <- sum(var_coexist)
    pred_coexist[c,4] <- sum(var_sp1)
    pred_coexist[c,5] <- sum(var_sp2)
    pred_coexist[c,6] <- sum(var_prior)
    
    c <- c+1
    
    coexist <- rbind(coexist, temp_coexist)
  }
}

#plots
#plot 1 - separate bacteria-bacteria and yeast-yeast competitions, superimposed on OG chesson plot
fun_default <- function(x) {x}
fun_inverse <- function(x) {1/x}
fun_vert <- function(x,y) {y <- 1}

colour_b_y <- c("no" = "#212121",
                "yes" = "#ff0000")

plot_chesson <- ggplot(data = coexist, mapping = aes(x = rho, y = f2_f1)) +
  theme_bw() +
  #geom_point(colour = '#212121', shape = 21, size = 1.6, stroke = 1.0, alpha = 0.5) +
  geom_point(aes(colour = V5), shape = 21, size = 1.6, stroke = 1.0, alpha = 0.5) +
  #scale_fill_manual(values = colour_b_y) +
  scale_color_manual(values = colour_b_y) +
  xlim(0, 2) + #truncating x axis at rho = 1, eliminates the area of priority effects as well as making the niche overlap close to 0 more obvious
  #ylim(0, 10) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary
plot_chesson
ggsave("~/Documents/figures_dump/chesson_black.pdf", width=6.4, height=6.4, units = "in")

#do a ttest on rho to see the difference between groups
t.test(subset(coexist, coexist$V5 == "yes")$rho,subset(coexist, coexist$V5 == "no")$rho, alternative = "t")
t.test(subset(coexist, coexist$V5 == "yes")$f2_f1,subset(coexist, coexist$V5 == "no")$f2_f1, alternative = "t")

#%coexisting - for b-b, b-y, y-y
b_b <- c(1,2,3,4,9,10,11,16,17,22)
b_y <- c(5,6,7,8,12,13,14,15,18,19,20,21,23,24,25,26,27,28,29,30)
y_y <- c(31,32,33,34,35,36)

pred_cox_tog <- data.frame(category=character(),
                           coexist = integer(),
                           S1_win = integer(),
                           S2_win = integer(),
                           priority = integer())

pred_cox_tog[1,1] <- 'b_b' 
pred_cox_tog[2,1] <- 'b_y' 
pred_cox_tog[3,1] <- 'y_y'
pred_cox_tog[4,1] <- 'all'
pred_cox_tog[1,2] <- sum(pred_coexist[b_b,]$coexist)/500
pred_cox_tog[1,3] <- sum(pred_coexist[b_b,]$S1_win)/500
pred_cox_tog[1,4] <- sum(pred_coexist[b_b,]$S2_win)/500
pred_cox_tog[1,5] <- sum(pred_coexist[b_b,]$priority)/500
pred_cox_tog[2,2] <- sum(pred_coexist[b_y,]$coexist)/1000
pred_cox_tog[2,3] <- sum(pred_coexist[b_y,]$S1_win)/1000
pred_cox_tog[2,4] <- sum(pred_coexist[b_y,]$S2_win)/1000
pred_cox_tog[2,5] <- sum(pred_coexist[b_y,]$priority)/1000
pred_cox_tog[3,2] <- sum(pred_coexist[y_y,]$coexist)/300
pred_cox_tog[3,3] <- sum(pred_coexist[y_y,]$S1_win)/300
pred_cox_tog[3,4] <- sum(pred_coexist[y_y,]$S2_win)/300
pred_cox_tog[3,5] <- sum(pred_coexist[y_y,]$priority)/300
pred_cox_tog[4,2] <- sum(pred_coexist$coexist)/1800
pred_cox_tog[4,3] <- sum(pred_coexist$S1_win)/1800
pred_cox_tog[4,4] <- sum(pred_coexist$S2_win)/1800
pred_cox_tog[4,5] <- sum(pred_coexist$priority)/1800

print(pred_cox_tog)

#make pie charts for the proportion of species that coexist (or not?)
category <- colnames(pred_coexist[,3:6])
colour_list <- c("#FDAE99","#882255","#DDCC77","#CC6677","#aa4499","#117733","#88ccee","#332288","#44aa99")
count <- 1
plot_list <- list()

for(i in 1:8){
  S1 = strain_names[i]
  for(j in 2:9){
    if(j <= i){
      plot_list[[length(plot_list)+1]] <- pNew
      plot_list[length(plot_list)] <- list(NULL)
      print(length(plot_list))
    }
    else{
      S2 = strain_names[j]
      cat(S1, ", ", S2, "\n")
      print(count)
      amount <- c(t(pred_coexist[count,3:6]))
      
      temp_data <- tibble(category = category,
                          amount = amount)
      
      pNew <- ggplot(temp_data, aes(x="", y=amount, fill=category)) +
        theme_classic() +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        #geom_text(aes(label = amount), position = position_stack(vjust=0.5), size = 2.5) +
        labs(x = NULL, y = NULL) +
        scale_fill_manual(values = c("coexist" = "#808080", 
                                     "S1_win" = colour_list[i], 
                                     "S2_win"= colour_list[j], 
                                     "priority" = "#202020")) +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position="none")
      plot_list[[length(plot_list)+1]] <- pNew
      count = count + 1
    }
  }
}

#plot all the bar plots in a grid
p_all <- plot_grid(
  plotlist = plot_list,
  ncol = 8,
  byrow = TRUE
)
p_all #1048X1048
ggsave("~/Documents/figures_dump/chesson_pie.pdf", width=10.48, height=10.48, units = "in")


#plot of all single and pairwise growth curves - supplementary figure
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177") #order of strains with strain IDs
colour_strainID <- c("163" = "#88CCEE",
                     "228" = "#332288",
                     "253" = "#117733",
                     "177" = "#44AA99",
                     "17B2" = "#FDAE99",
                     "0092a" = "#882255",
                     "232" = "#DDCC77",
                     "460" = "#CC6677",
                     "550" = "#AA4499") #colour names with strain IDs
 
#plot 1: single species growth curves
GRdata_plot <- read.csv("/Users/kasturilele/Documents/community/GCobs_single_normalized.csv", header = T)
GRdata_plot$Strain <- factor(GRdata_plot$Strain, levels=strain_order)
GRdata_model_single <- read.csv("/Users/kasturilele/Documents/cluster/allSingSpecTraj_5-3.csv", header = T)
colnames(GRdata_model_single) <- c("Strain", "Replicate", "Time", "CFU_strain_1")
GRdata_model_single <- subset(GRdata_model_single, Replicate < 20)

#single species facet plot - two versions for final figure (horizontal and vertical)
#horizontal - 14.75X1.55, vertical - 2.25X11.05
ggplot() + 
  theme_bw()+
  geom_point(data = GRdata_plot, mapping = aes(x = Time,y = CFU_Strain_1, group = Replicate, colour = Strain)) + #plots the points for strain 1
  geom_line(data = GRdata_model_single, mapping = aes(x = Time, y = CFU_strain_1, group = Replicate, colour = Strain), alpha = 0.25) +
  scale_color_manual(values = colour_strainID) +
  ylim(0,4e07)+
  facet_wrap(vars(factor(Strain, levels = strain_order)), nrow =1)+
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none",
        strip.text.x = element_blank()) 
ggsave("~/Documents/for_paper/figures_dump/1_single_horz.pdf", width=14.75, height=1.55, units = "in")

ggplot() + 
  theme_bw()+
  geom_point(data = GRdata_plot, mapping = aes(x = Time,y = CFU_Strain_1, group = Replicate, colour = Strain)) + #plots the points for strain 1
  geom_line(data = GRdata_model_single, mapping = aes(x = Time, y = CFU_strain_1, group = Replicate, colour = Strain), alpha = 0.25) +
  scale_color_manual(values = colour_strainID) +
  ylim(0,4e07)+
  facet_wrap(vars(factor(Strain, levels = strain_order)), nrow =9)+
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none",
        strip.text.x = element_blank()) 
ggsave("~/Documents/for_paper/figures_dump/1_single_vert.pdf", width=2.25, height=11.05, units = "in")

# plot 2: pairwise growth curves 
GRdatapw <-read.csv("/Users/kasturilele/Documents/community/GCobs_paired_normalized.csv", header = T)
GRdatamodel <- read.csv("/Users/kasturilele/Documents/cluster/allPairSpecTraj_6-3.csv", header = T)
GRdata_model_paired = data.frame(matrix(nrow = 0, ncol = 6))
colnames(GRdata_model_paired) <- c("Strain_1","Strain_2","Replicate","Time","CFU_strain_1","CFU_strain_2")

GRdatapw$Strain_1 <- factor(GRdatapw$Strain_1, levels=strain_order)
GRdatapw$Strain_2 <- factor(GRdatapw$Strain_2, levels=strain_order)

i <- 1
for(i in 1:8){
  for(j in (i+1):9){
    temp_paired = subset(GRdatamodel, Strain_1 == strain_order[i]&Strain_2 == strain_order[j])
    if(nrow(temp_paired) == 0){
      temp_paired = subset(GRdatamodel, Strain_2 == strain_order[i]&Strain_1 == strain_order[j])
      colnames(temp_paired) <- c("Strain_2","Strain_1","Replicate","Time","CFU_strain_2","CFU_strain_1")
    }
    GRdata_model_paired <- rbind(GRdata_model_paired,temp_paired)
  }
}
GRdata_model_paired$Strain_1 <- factor(GRdata_model_paired$Strain_1, levels=strain_order)
GRdata_model_paired$Strain_2 <- factor(GRdata_model_paired$Strain_2, levels=strain_order)
GRdata_model_paired <- subset(GRdata_model_paired, Replicate < 20)

ggplot(data = GRdatapw, mapping = aes(x = time, group = Replicate)) +
  theme_bw()+
  geom_point(aes(y = CFU_Strain_1, colour = Strain_1)) + #plots the points for strain 1
  geom_line(data = GRdata_model_paired, mapping = aes(x = Time, y = CFU_strain_1, group = Replicate, colour = Strain_1), alpha = 0.25) +
  geom_point(aes(y = CFU_Strain_2, colour = Strain_2)) + #plots strain 2
  geom_line(data = GRdata_model_paired, mapping = aes(x = Time, y = CFU_strain_2, group = Replicate, colour = Strain_2), alpha = 0.25) +
  scale_color_manual(values = colour_strainID) +
  #scale_y_log10() + #sets the y scale to a log scale instead of a linear scale
  facet_grid(Strain_2~Strain_1) +
  #, scale = "free_y") + #plots all the strains separately
  labs(x = "time (hrs)", y = "total CFU/well") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")

ggsave("~/Documents/for_paper/figures_dump/1_pairwise.pdf", width=13.44, height=10.11, units = "in")


#--------- figure 3: multi-species model prediction ---------
library(ggplot2)
library(cowplot)
library(grid)
library(dplyr)
library(vegan)
library(indicspecies)
library(tidyr)
#open the file 
new_data_main <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_OG_model_288.csv", header = T, row.names = 1))
new_data_plot <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_OG_model_288.csv", header = T))
new_data_main <- round(new_data_main)
new_data_plot[,c(2:10)] <- round(new_data_plot[,c(2:10)]) #to get rid of the one spurious negative value that somehow snuck into the data

new_data_main[new_data_main < 100000] = 0 #remove species below detection limit

colour_names <- c("F.sanfranciscensis" = "#FDAE99",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99")

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum",
                  "C.paralimentarius","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

#processing the data so that it is relative abundance instead of absolute abundance
new_data <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(new_data) <- colnames(new_data_main)
for(i in 1:100){
  temp_data <- new_data_main[i,]
  temp_data <- temp_data/sum(temp_data)
  new_data <- rbind(new_data,temp_data)
}

#new community assembly data - 
compare_df <- as.data.frame(read.csv("~/Documents/cluster/new_CFU_relative.csv", header = T, row.names = 1))
compare_long <- as.data.frame(read.csv("~/Documents/cluster/new_CFU.csv", header = T))

compare_df <- compare_df[1:10,]
compare_long <- compare_long[1:10,c(1,3:11)]

Sample_ID <- c(0:9)
Source <- compare_df$Origin

compare_meta <- as.data.frame(cbind(Sample_ID, Source))
compare_df <- compare_df[,2:10]

Sample_ID <- new_data_plot$Population
Source <- rep('M', times = 100)
model_meta <- as.data.frame(cbind(Sample_ID,Source))

both_df <- rbind(new_data,compare_df)
both_meta <- rbind(model_meta,compare_meta)

#community composition
#distance matrix - bray-curtis dissimilarity
distance_model <- vegdist(new_data, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_model <- hclust(distance_model, method = "ward.D2")
plot(hc_model)

dd_model <- as.dendrogram(hc_model)
pdf("~/Documents/for_paper/figures_dump/2_model_dend.pdf", width=13.6, height=3.20) # The height of the plot in inches
plot(dd_model) #1242X521
dev.off()

#make a data frame with the cluster order
model_clustorder<-hc_model$labels[c(hc_model$order)] 
model_clustorder<-as.data.frame(model_clustorder)
model_clustorder<-cbind(model_clustorder,rownames(model_clustorder))
colnames(model_clustorder) <- c("Population","order")
model_clusters <- model_clustorder

# change data from wide format to long format
colIDs <- colnames(new_data)
new_data_long <- new_data_plot %>% pivot_longer(cols=all_of(colIDs),
                                                names_to='Species',
                                                values_to='Abundance')
new_data_long <- merge(new_data_long, model_clusters, by = 'Population')

new_data_long$order<-as.numeric(new_data_long$order)
new_data_long<-arrange(new_data_long, order)

p_ordered <- ggplot(data = new_data_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = factor(Species, levels = strain_names))) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_ordered #1242X261
ggsave("~/Documents/for_paper/figures_dump/2_model_clustered.pdf", width=12.42, height=2.60, units = "in")

#make a data frame with average model and observed community
avg_both <- data.frame(Species=character(),
                       avg_model=double(),
                       avg_obs=double())

#measuring average proportion of each species
i <- 1
for(i in 1:9){
  species <- new_data[,i]
  avg_species <- mean(species)
  print(colnames(new_data)[i])
  print(avg_species*100)
  avg_both[i,1] <- colnames(new_data)[i]
  avg_both[i,2] <- avg_species
}

#cluster for obs data
distance_obs <- vegdist(compare_df, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_obs <- hclust(distance_obs, method = "ward.D2")

dd_obs <- as.dendrogram(hc_obs)

pdf("~/Documents/for_paper/figures_dump/2_obs_dend.pdf", width=2.42, height=2.60) # The height of the plot in inches
plot(dd_obs) #1242X521
dev.off()

#make a data frame with the cluster order
obs_clustorder<-hc_obs$labels[c(hc_obs$order)] 
obs_clustorder<-as.data.frame(obs_clustorder)
obs_clustorder<-cbind(obs_clustorder,rownames(obs_clustorder))
colnames(obs_clustorder) <- c("Population","order")
obs_clusters <- obs_clustorder

#plot

# change data from wide format to long format
colIDs <- colnames(compare_df)
compare_long <- compare_long %>% pivot_longer(cols=all_of(colIDs),
                                              names_to='Species',
                                              values_to='Abundance')
compare_long <- merge(compare_long, obs_clusters, by = 'Population')

compare_long$order<-as.numeric(compare_long$order)
compare_long<-arrange(compare_long, order)

p_obs <- ggplot(data = compare_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = factor(Species, levels = strain_names))) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_obs #1242X261
ggsave("~/Documents/for_paper/figures_dump/observed_all.pdf", width=2.42, height=2.60, units = "in")

#measuring average proportion of each species
i <- 1
for(i in 1:9){
  species <- compare_df[,i]
  avg_species <- mean(species)
  print(colnames(compare_df)[i])
  print(avg_species*100)
  avg_both[i,3] <- avg_species#avg_both was made earlier for the model avgs
}

#shannon,bray-curtis,NMDS for the model+obs combined

#1 - diversity indices - alpha diversity
both_meta$Shannon <- diversity(both_df,"shannon")

#plot data wrt treatment (model/observed). do a bunch of data shifting to merge with model clusters in order to colour model data
model_shannon_temp <- both_meta[1:100,]
colnames(model_shannon_temp) <- c("Population","Source","Shannon")
model_shannon <- merge(model_shannon_temp, model_clusters, by = 'Population')
obs_shannon <- both_meta[101:110,]
colnames(obs_shannon) <- c("Population","Source","Shannon")
obs_shannon$order <- c(0,0,0,0,0,0,0,0,0,0)
shannon <- rbind(model_shannon, obs_shannon)

both_shannon <- ggplot(shannon, aes(x = Source, y = Shannon)) +
  theme_bw() +
  geom_boxplot(color = '#888888', outlier.shape = NA) +
  geom_jitter(data = model_shannon, shape = 16, size = 5.0, alpha = 0.60, aes(colour = order), width = 0.25, height = 0) +
  scale_color_viridis_d(option = "B") + 
  geom_jitter(data = obs_shannon, shape = 17, size = 3.0, colour = "#000000", width = 0.2, height = 0) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
both_shannon
ggsave("~/Documents/for_paper/figures_dump/2_both_shannon.pdf", width=5, height=7, units = "in")
t.test(Shannon ~ Source, both_meta) #t-test for finding out differences between shannon indices

#2 - community composition

#distance matrix - bray-curtis dissimilarity
both_distance <- as.matrix(vegdist(both_df, "bray")) #distance matrix
#do a permanova
both_adonis <- adonis2(both_distance~Source, data = both_meta, permutations = 1000)
both_adonis

#do nmds to visualize this above stuff
nmds_both <- metaMDS(both_distance, k = 2, maxit = 999, try = 500)
plot(nmds_both)
#save the coordinates
coords_both <- as.data.frame(nmds_both$points)

#add columns to data frame 
coords_both$Source = both_meta$Source

#separate model/obs data to create graphs coloured by hierarchical clustering of model
coords_nmds_model <- coords_both[1:100,]
coords_nmds_obs <- coords_both[101:110,]
coords_nmds_model$Population <- rownames(coords_nmds_model)
nmds_plot_data <- merge(coords_nmds_model, model_clusters, by = 'Population')
nmds_plot_data$order <- as.numeric(nmds_plot_data$order)/100

#plot nmds coloured by gradient
p_grad = ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 16, size = 5.0, alpha = 0.60, aes(colour = order)) +
  geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 17, size = 3.0, colour = "#000000") +
  scale_color_viridis_c(option = "B") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_grad

#overlaying species scores (wascores)
sp <- as.data.frame(wascores(x = nmds_both$points, w = both_df, expand = TRUE))
p_grad_sp <- p_grad +
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = sp, size =1, alpha = 0.5, colour = "#4d4d4d") +
  geom_text(data = sp, aes(x = MDS1, y = MDS2), alpha = 0.75,
            label = row.names(sp))
p_grad_sp


ggsave("~/Documents/for_paper/figures_dump/2_nmds_both_spscores.pdf", width=6.40, height=6.40, units = "in")

#3 - model vs experiment ranked by abundance
#3a model data ranked by abundance
ranked_model <- data.frame(new_data) %>%
  rowwise() %>% 
  do(data.frame(t(rank(-unlist(.)))))

rownames(ranked_model) <- rownames(new_data)

#sum how many communities for each species have which rank
sum_ranked_model <- data.frame(Species=character(),
                               rank=double(),
                               sum_rank=integer())

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
counter <- 1
i <- 1
for(i in 1:9){
  ranked_col <- ranked_model[[i]]
  ranks <- unique(ranked_col)
  L <- length(ranks)
  j <- 1
  for(j in 1:L){
    rank_num <- which(ranked_col == ranks[j])
    sum_ranked_model[counter,1] <- strain_names[i]
    sum_ranked_model[counter,2] <- ranks[j]
    sum_ranked_model[counter,3] <- length(rank_num)
    counter <- counter + 1
  }
}

#3b observed data ranked by abundance
ranked_obs <- data.frame(compare_df) %>%
  rowwise() %>% 
  do(data.frame(t(rank(-unlist(.)))))

rownames(ranked_obs) <- rownames(compare_df)

#sum how many communities for aeach species have which rank
sum_ranked_obs <- data.frame(Species=character(),
                             rank=double(),
                             sum_rank=integer())

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
counter <- 1
i <- 1
for(i in 1:9){
  ranked_col <- ranked_obs[[i]]
  ranks <- unique(ranked_col)
  L <- length(ranks)
  j <- 1
  for(j in 1:L){
    rank_num <- which(ranked_col == ranks[j])
    sum_ranked_obs[counter,1] <- strain_names[i]
    sum_ranked_obs[counter,2] <- ranks[j]
    sum_ranked_obs[counter,3] <- length(rank_num)
    counter <- counter + 1
  }
}

#3c plotting model and observed ranks against each other
sum_ranked_model$sum_rank = sum_ranked_model$sum_rank/100
Source <- rep('M', nrow(sum_ranked_model))
sum_ranked_model <- cbind(sum_ranked_model,Source)
sum_ranked_obs$sum_rank = sum_ranked_obs$sum_rank/10
Source <- rep('O', nrow(sum_ranked_obs))
sum_ranked_obs <- cbind(sum_ranked_obs,Source)
sum_ranked <- rbind(sum_ranked_model,sum_ranked_obs)
sum_ranked$Species <- factor(sum_ranked$Species, levels = strain_names)
p3 <- ggplot(data = sum_ranked, mapping = aes(x = Source, y = rank, colour = Species, size = sum_rank))+
  theme_bw() +
  geom_point() +
  scale_colour_manual(values = colour_names) +
  scale_y_reverse() +
  facet_wrap(~ Species, nrow = 1) +
  theme(#panel.grid.minor.x = element_blank(),
    #axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none")
p3
ggsave("~/Documents/for_paper/figures_dump/2_bubble_both.pdf", width=10.24, height=3.2, units = "in")

#3d doing a spearman correlation on average model and observed rank
spear_both <- cor.test(avg_both$avg_model, avg_both$avg_obs, method = 'spearm', alternative = 'two.sided', exact = FALSE)
spear_both

#indicator species analysis
indval_both <- multipatt(both_df, both_meta$Source,
                         control = how(nperm=999))
summary(indval_both, alpha = 1)

#plot of indicator species - supplementary figure
relAbdPlot <- list()
i <- 1
colours <- c( "#FDAE99","#882255","#DDCC77","#aa4499","#CC6677","#117733","#88ccee","#332288","#44aa99")
for(i in 1:9){
  coords_sp <- cbind(coords_both, both_df[,i])
  colnames(coords_sp) <- c('MDS1','MDS2','Source','rel_abd')
  coords_sp_plot <- coords_sp[which(coords_sp$rel_abd > 0),]
  cols <- colours[i]
  #plot nmds coloured by gradient
  p_sp = ggplot(coords_sp_plot, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_point(shape = 21, stroke = 1.0, colour = cols, aes(size = rel_abd)) +
    geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 1.0, colour = "#000000", alpha = 0.5) +
    geom_point(data = coords_nmds_model, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 1.0, colour = "#888888", alpha = 0.25) +
    coord_cartesian(xlim = c(-0.3,0.25), ylim = c(-0.2, 0.25)) +
    #scale_color_manual(values = colour_names) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')
  p_sp
  relAbdPlot[[length(relAbdPlot)+1]] <- p_sp
}

relAbdAll <- plot_grid(
  plotlist = relAbdPlot,
  ncol = 3,
  byrow = TRUE
)
relAbdAll

ggsave("~/Documents/for_paper/figures_dump/2_nmds_relsp.pdf", width=6.4, height=6.4, units = "in")

#--------- figure 6: leave-one-out data ---------

library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(cowplot)

colour_names <- c("F.sanfranciscensis" = "#FDAE99",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99")

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

#comparison between all observed leave-one-outs
plot_list <- list()
obs_data <- as.data.frame(read.csv("/Users/kasturilele/Documents/cluster/new_CFU_relative.csv", header = T))
obs_data <- obs_data[rowSums(obs_data[,3:11])>0,]

#use this loop to generate and save all the relative abundance and hierarchical clustering plots for communities
i <- 0
for(i in 0:8){
  #open each new file
  filename <- paste("output_", i,".csv",sep = "" )
  filesave_dend <- paste("~/Documents/figures_dump/LOO/dend_",i, ".pdf", sep = "")
  filesave_clust <- paste("~/Documents/figures_dump/LOO/clustered_",i, ".pdf", sep = "")
  #open the data
  LOO_plot <- subset(obs_data, Origin == strain_names[i+1])[,c(1,3:11)]
  LOO_data <- LOO_plot[,2:10]
  rownames(LOO_data) <- LOO_plot$Population
  
  dm_model_data <- vegdist(LOO_data, "bray")
  hc_model <- hclust(dm_model_data, method = "ward.D2")
  
  dd_model <- as.dendrogram(hc_model)
  pdf(filesave_dend, width=2.42, height=2.60) # The height of the plot in inches
  plot(dd_model)
  dev.off()
  
  #make a data frame with the cluster order
  model_clustorder<-hc_model$labels[c(hc_model$order)] 
  model_clustorder<-as.data.frame(model_clustorder)
  model_clustorder<-cbind(model_clustorder,rownames(model_clustorder))
  colnames(model_clustorder) <- c("Population","order")
  model_clusters <- model_clustorder
  
  #plot
  
  # change data from wide format to long format
  colIDs <- colnames(LOO_data)
  LOO_long <- LOO_plot %>% pivot_longer(cols=all_of(colIDs),
                                        names_to='Species',
                                        values_to='Abundance')
  
  LOO_long <- merge(LOO_long, model_clusters, by = 'Population')
  
  LOO_long$order<-as.numeric(LOO_long$order)
  LOO_long<-arrange(LOO_long, order)
  
  p_ordered <- ggplot(data = LOO_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = factor(Species, levels = strain_names))) +
    theme_minimal() +
    geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
    scale_fill_manual(values = colour_names) +
    labs(title = strain_names[i+1]) + 
    theme(axis.text = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, size = 16),
          legend.position = 'none')
  
  p_ordered
  ggsave(filesave_clust, width=2.42, height=2.60, units = "in")
  
}

#-----------figure 6: comparison between all observed communities (9-species and LOO)-----------
Sample_ID <- obs_data$Population
Source <- obs_data$Origin
obs_df <- obs_data[,3:11]
rownames(obs_df) <- obs_data$Population
obs_meta <- as.data.frame(cbind(Sample_ID,Source))

#analysis of observed community composition

#distance matrix - bray-curtis dissimilarity
obs_distance <- as.matrix(vegdist(obs_df, "bray")) #distance matrix
#do a permanova
obs_adonis <- adonis2(obs_distance~Source, data = obs_meta, permutations = 1000)
obs_adonis

#do nmds to visualize this above stuff
nmds_obs <- metaMDS(obs_distance, k = 2, maxit = 999, try = 500)
plot(nmds_obs)
#save the coordinates
coords_obs <- as.data.frame(nmds_obs$points)

#add columns to data frame 
coords_obs$Source = obs_meta$Source

#separate model/obs data to create graphs coloured by hierarchical clustering of model

#plot nmds coloured by gradient
p_nmds = ggplot(coords_obs, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 16, size = 5.0, alpha = 0.60, aes(colour = Source)) +
  scale_color_manual(values = colour_names) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_nmds

#messing around with species scores (wascores vs envdist)
sp <- as.data.frame(wascores(x = nmds_obs$points, w = obs_df, expand = TRUE))
p_nmds_sp <- p_nmds +
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = sp, size =1, alpha = 0.5, colour = "#4d4d4d") +
  geom_text(data = sp, aes(x = MDS1, y = MDS2), alpha = 0.75,
            label = row.names(sp))
p_nmds_sp
ggsave("~/Documents/figures_dump/LOO/nmds_all_obs.pdf", width=6.40, height=6.40, units = "in")

# mindist comparisons for model:obs for all vs LOOs
setwd("/Users/kasturilele/Documents/community")
mindist_combi <- data.frame("all" = double(),
                            "F.sanfranciscensis" = double(),
                            "L.brevis" = double(),
                            "L.plantarum" = double(),
                            "A.malorum" = double(),
                            "C.paralimentarius" = double(),
                            "S.cerevisiae" = double(),
                            "W.anomalus" = double(),
                            "K.humilis" = double(),
                            "K.servazzii" = double())

for(n in 0:8){
  #open each new file
  filename <- paste("output_ST_", n,".csv",sep = "" )
  #open the data
  LOO_model_main <- read.csv(filename, header = T, row.names = 1)
  
  #processing the data so that it is relative abundance instead of absolute abundance
  LOO_model <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(LOO_model) <- colnames(LOO_model_main)
  for(i in 1:100){
    temp_data <- LOO_model_main[i,]
    temp_data <- temp_data/sum(temp_data)
    LOO_model <- rbind(LOO_model,temp_data)
  }
  
  LOO_obs <- subset(obs_data, Origin == strain_names[n+1])[,3:11]
  
  Sample_ID <- rownames(LOO_obs)
  Source <- subset(obs_data, Origin == strain_names[n+1])$Origin
  
  compare_meta <- as.data.frame(cbind(Sample_ID, Source))
  
  Sample_ID <- rownames(LOO_model)
  Source <- rep('M', times = 100)
  model_meta <- as.data.frame(cbind(Sample_ID,Source))
  
  both_df <- rbind(LOO_model,LOO_obs)
  both_meta <- rbind(model_meta,compare_meta)
  
  both_distance <- as.matrix(vegdist(both_df, "bray")) #distance matrix
  
  both_distance_sub <- as.data.frame(both_distance[101:nrow(both_df),1:100])
  ModelID <- colnames(both_distance_sub)
  
  i = 1
  for(i in 1:100){
    vector_temp <- both_distance_sub[i]
    minval <- min(vector_temp)
    mindist_combi[i,n+2] <- minval
  }
}
#adding the results for 9-species community
LOO_model_main <- read.csv("/Users/kasturilele/Documents/cluster/endpoints_new_serial_model.csv", header = T, row.names = 1)

#processing the data so that it is relative abundance instead of absolute abundance
LOO_model <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(LOO_model) <- colnames(LOO_model_main)
for(i in 1:100){
  temp_data <- LOO_model_main[i,]
  temp_data <- temp_data/sum(temp_data)
  LOO_model <- rbind(LOO_model,temp_data)
}

LOO_obs <- subset(obs_data, Origin == "all")[,3:11]

Sample_ID <- rownames(LOO_obs)
Source <- subset(obs_data, Origin == "all")$Origin

compare_meta <- as.data.frame(cbind(Sample_ID, Source))

Sample_ID <- rownames(LOO_model)
Source <- rep('M', times = 100)
model_meta <- as.data.frame(cbind(Sample_ID,Source))

both_df <- rbind(LOO_model,LOO_obs)
both_meta <- rbind(model_meta,compare_meta)

both_distance <- as.matrix(vegdist(both_df, "bray")) #distance matrix

both_distance_sub <- as.data.frame(both_distance[101:nrow(both_df),1:100])
ModelID <- colnames(both_distance_sub)

i = 1
for(i in 1:100){
  vector_temp <- both_distance_sub[i]
  minval <- min(vector_temp)
  mindist_combi[i,1] <- minval
}
mindist_long <- mindist_combi %>% pivot_longer(everything(),
                                               names_to='Species',
                                               values_to='mindist')
#set order of factors
source_order <- unique(mindist_long$Species)
mindist_long$Species <- factor(mindist_long$Species, levels = source_order)

#histogram version 7.2X2.4
plot_mindist <- ggplot(data = mindist_combi) +
  theme_bw() +
  geom_density(aes(all),colour = "#757575", adjust = 0.75, linewidth = 3) +
  geom_density(aes(F.sanfranciscensis),colour = "#FDAE99", adjust = 0.75) +
  geom_density(aes(L.brevis),colour = "#882255", adjust = 0.75) +
  geom_density(aes(L.plantarum),colour = "#DDCC77", adjust = 0.75) +
  geom_density(aes(A.malorum),colour = "#CC6677", adjust = 0.75) +
  geom_density(aes(C.paralimentarius),colour = "#AA4499", adjust = 0.75) +
  geom_density(aes(S.cerevisiae),colour = "#117733", adjust = 0.75) +
  geom_density(aes(W.anomalus),colour = "#88ccee", adjust = 0.75) +
  geom_density(aes(K.humilis),colour = "#332288", adjust = 0.75) +
  geom_density(aes(K.servazzii),colour = "#44aa99", adjust = 0.75) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
plot_mindist
ggsave("~/Documents/figures_dump/LOO/all_mindist_serial.pdf", width=10.24, height=3.60, units = "in")

#------------- figure 4: serial bottleneck model ------------
library(ggplot2)
library(cowplot)
library(grid)
library(dplyr)
library(vegan)
library(indicspecies)
library(tidyr)

colour_names <- c("F.sanfranciscensis" = "#FDAE99",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99")

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum",
                  "C.paralimentarius","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

#open the original model file 
OG_data_main <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_OG_model_288.csv", header = T, row.names = 1))
OG_data_plot <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_OG_model_288.csv", header = T))
OG_data_main <- round(OG_data_main)
OG_data_plot[,c(2:10)] <- round(OG_data_plot[,c(2:10)]) #to get rid of the one spurious negative value that somehow snuck into the data

OG_data_main[OG_data_main < 100000] = 0 #remove species below detection limit

#processing the data so that it is relative abundance instead of absolute abundance
OG_data <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(OG_data) <- colnames(OG_data_main)
for(i in 1:100){
  temp_data <- OG_data_main[i,]
  temp_data <- temp_data/sum(temp_data)
  OG_data <- rbind(OG_data,temp_data)
}

#open the file with the serial dilution model results
new_data_main <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_serial_model.csv", header = T, row.names = 1))
new_data_plot <- as.data.frame(read.csv("~/Documents/cluster/endpoints_new_serial_model.csv", header = T))
new_data_main <- round(new_data_main)
new_data_plot[,c(2:10)] <- round(new_data_plot[,c(2:10)]) #to get rid of the one spurious negative value that somehow snuck into the data

new_data_main[new_data_main < 100000] = 0 #remove species below detection limit

#processing the data so that it is relative abundance instead of absolute abundance
new_data <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(new_data) <- colnames(new_data_main)
for(i in 1:100){
  temp_data <- new_data_main[i,]
  temp_data <- temp_data/sum(temp_data)
  new_data <- rbind(new_data,temp_data)
}

#hierarchical clustering of serial transfer data
#distance matrix - bray-curtis dissimilarity
distance_model <- vegdist(new_data, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_model <- hclust(distance_model, method = "ward.D2")
plot(hc_model)

dd_model <- as.dendrogram(hc_model)
pdf("~/Documents/for_paper/figures_dump/4_serial_dend.pdf", width=13.6, height=3.20) # The height of the plot in inches
plot(dd_model) #1242X521
dev.off()

#make a data frame with the cluster order
model_clustorder<-hc_model$labels[c(hc_model$order)] 
model_clustorder<-as.data.frame(model_clustorder)
model_clustorder<-cbind(model_clustorder,rownames(model_clustorder))
colnames(model_clustorder) <- c("Population","order")
model_clusters <- model_clustorder

# change data from wide format to long format
colIDs <- colnames(new_data)
new_data_long <- new_data_plot %>% pivot_longer(cols=all_of(colIDs),
                                                names_to='Species',
                                                values_to='Abundance')
new_data_long <- merge(new_data_long, model_clusters, by = 'Population')

new_data_long$order<-as.numeric(new_data_long$order)
new_data_long<-arrange(new_data_long, order)

p_ordered <- ggplot(data = new_data_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = factor(Species, levels = strain_names))) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_ordered #1242X261
ggsave("~/Documents/for_paper/figures_dump/4_serial_clustered.pdf", width=12.42, height=2.60, units = "in")

#additional analyses on serial bottleneck model vs experiment
compare_df <- as.data.frame(read.csv("~/Documents/cluster/new_CFU_relative.csv", header = T, row.names = 1))
compare_df <- compare_df[1:10,]


Sample_ID <- c(0:9)
Source <- compare_df$Origin

compare_meta <- as.data.frame(cbind(Sample_ID, Source))
compare_df <- compare_df[,2:10]

Sample_ID <- OG_data_plot$Population
Source <- rep('M', times = 100)
model_meta <- as.data.frame(cbind(Sample_ID,Source))

both_df_OG <- rbind(OG_data,compare_df)
both_df_new <- rbind(new_data,compare_df)
both_meta <- rbind(model_meta,compare_meta)

#making mindist plots
mindist_both <- data.frame("OG" = double(),
                           "new" = double())

both_distance <- as.matrix(vegdist(both_df_OG, "bray")) #distance matrix

both_distance_sub <- as.data.frame(both_distance[101:110,1:100])

i = 1
for(i in 1:100){
  vector_temp <- both_distance_sub[i]
  minval <- min(vector_temp)
  mindist_both[i,1] <- minval
}

both_distance <- as.matrix(vegdist(both_df_new, "bray")) #distance matrix

both_distance_sub <- as.data.frame(both_distance[101:110,1:100])

i = 1
for(i in 1:100){
  vector_temp <- both_distance_sub[i]
  minval <- min(vector_temp)
  mindist_both[i,2] <- minval
}

mindist_long <- mindist_both %>% pivot_longer(everything(),
                                              names_to='origin',
                                              values_to='mindist')

#histogram version 7.2X2.4
plot_mindist <- ggplot(data = mindist_both) +
  theme_bw() +
  geom_density(aes(OG),colour = "#cccccc", adjust = 0.75) +
  geom_density(aes(new),colour = "#111111", adjust = 0.75) +
  xlim(0,0.5) +
  ylim(0,30) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
plot_mindist
ggsave("~/Documents/for_paper/revisions/figdump/mindist.pdf", width=7.20, height=2.40, units = "in")
#redoing analyses with serial transfer data
#shannon,bray-curtis,NMDS for the model+obs combined

#1 - diversity indices - alpha diversity
both_meta$Shannon <- diversity(both_df_new,"shannon")

#plot data wrt treatment (model/observed). do a bunch of data shifting to merge with model clusters in order to colour model data
model_shannon_temp <- both_meta[1:100,]
colnames(model_shannon_temp) <- c("Population","Source","Shannon")
model_shannon <- merge(model_shannon_temp, model_clusters, by = 'Population')
obs_shannon <- both_meta[101:110,]
colnames(obs_shannon) <- c("Population","Source","Shannon")
obs_shannon$order <- c(0,0,0,0,0,0,0,0,0,0)
shannon <- rbind(model_shannon, obs_shannon)

both_shannon <- ggplot(shannon, aes(x = Source, y = Shannon)) +
  theme_bw() +
  geom_boxplot(color = '#888888', outlier.shape = NA) +
  geom_jitter(data = model_shannon, shape = 16, size = 5.0, alpha = 0.60, aes(colour = order), width = 0.25, height = 0) +
  scale_color_viridis_d(option = "B") + 
  geom_jitter(data = obs_shannon, shape = 17, size = 3.0, colour = "#000000", width = 0.2, height = 0) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
both_shannon
ggsave("~/Documents/for_paper/figures_dump/4_serial_shannon.pdf", width=5, height=7, units = "in")
t.test(Shannon ~ Source, both_meta) #t-test for finding out differences between shannon indices

#2 - community composition

#distance matrix - bray-curtis dissimilarity
both_distance <- as.matrix(vegdist(both_df_new, "bray")) #distance matrix
#do a permanova
both_adonis <- adonis2(both_distance~Source, data = both_meta, permutations = 1000)
both_adonis

#do nmds to visualize this above stuff
nmds_both <- metaMDS(both_distance, k = 2, maxit = 999, try = 500)
plot(nmds_both)
#save the coordinates
coords_both <- as.data.frame(nmds_both$points)

#add columns to data frame 
coords_both$Source = both_meta$Source

#separate model/obs data to create graphs coloured by hierarchical clustering of model
coords_nmds_model <- coords_both[1:100,]
coords_nmds_obs <- coords_both[101:110,]
coords_nmds_model$Population <- rownames(coords_nmds_model)
nmds_plot_data <- merge(coords_nmds_model, model_clusters, by = 'Population')
nmds_plot_data$order <- as.numeric(nmds_plot_data$order)/100

#plot nmds coloured by gradient
p_grad = ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 16, size = 5.0, alpha = 0.60, aes(colour = order)) +
  geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 17, size = 3.0, colour = "#000000") +
  scale_color_viridis_c(option = "B") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_grad

#overlaying species scores (wascores)
sp <- as.data.frame(wascores(x = nmds_both$points, w = both_df_new, expand = TRUE))
p_grad_sp <- p_grad +
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = sp, size =1, alpha = 0.5, colour = "#4d4d4d") +
  geom_text(data = sp, aes(x = MDS1, y = MDS2), alpha = 0.75,
            label = row.names(sp))
p_grad_sp


ggsave("~/Documents/for_paper/figures_dump/4_nmds_serial_spscores.pdf", width=6.40, height=6.40, units = "in")

#3 - model vs experiment ranked by abundance
#3a model data ranked by abundance
ranked_model <- data.frame(new_data) %>%
  rowwise() %>% 
  do(data.frame(t(rank(-unlist(.)))))

rownames(ranked_model) <- rownames(new_data)

#sum how many communities for each species have which rank
sum_ranked_model <- data.frame(Species=character(),
                               rank=double(),
                               sum_rank=integer())

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
counter <- 1
i <- 1
for(i in 1:9){
  ranked_col <- ranked_model[[i]]
  ranks <- unique(ranked_col)
  L <- length(ranks)
  j <- 1
  for(j in 1:L){
    rank_num <- which(ranked_col == ranks[j])
    sum_ranked_model[counter,1] <- strain_names[i]
    sum_ranked_model[counter,2] <- ranks[j]
    sum_ranked_model[counter,3] <- length(rank_num)
    counter <- counter + 1
  }
}

#3b observed data ranked by abundance
ranked_obs <- data.frame(compare_df) %>%
  rowwise() %>% 
  do(data.frame(t(rank(-unlist(.)))))

rownames(ranked_obs) <- rownames(compare_df)

#sum how many communities for aeach species have which rank
sum_ranked_obs <- data.frame(Species=character(),
                             rank=double(),
                             sum_rank=integer())

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
counter <- 1
i <- 1
for(i in 1:9){
  ranked_col <- ranked_obs[[i]]
  ranks <- unique(ranked_col)
  L <- length(ranks)
  j <- 1
  for(j in 1:L){
    rank_num <- which(ranked_col == ranks[j])
    sum_ranked_obs[counter,1] <- strain_names[i]
    sum_ranked_obs[counter,2] <- ranks[j]
    sum_ranked_obs[counter,3] <- length(rank_num)
    counter <- counter + 1
  }
}

#3c plotting model and observed ranks against each other
sum_ranked_model$sum_rank = sum_ranked_model$sum_rank/100
Source <- rep('M', nrow(sum_ranked_model))
sum_ranked_model <- cbind(sum_ranked_model,Source)
sum_ranked_obs$sum_rank = sum_ranked_obs$sum_rank/10
Source <- rep('O', nrow(sum_ranked_obs))
sum_ranked_obs <- cbind(sum_ranked_obs,Source)
sum_ranked <- rbind(sum_ranked_model,sum_ranked_obs)
sum_ranked$Species <- factor(sum_ranked$Species, levels = strain_names)
p3 <- ggplot(data = sum_ranked, mapping = aes(x = Source, y = rank, colour = Species, size = sum_rank))+
  theme_bw() +
  geom_point() +
  scale_colour_manual(values = colour_names) +
  scale_y_reverse() +
  facet_wrap(~ Species, nrow = 1) +
  theme(#panel.grid.minor.x = element_blank(),
    #axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none")
p3
ggsave("~/Documents/for_paper/figures_dump/4_bubble_both.pdf", width=10.24, height=3.2, units = "in")

#3d doing a spearman correlation on average model and observed rank
spear_both <- cor.test(avg_both$avg_model, avg_both$avg_obs, method = 'spearm', alternative = 'two.sided', exact = FALSE)
spear_both

#indicator species analysis
indval_both <- multipatt(both_df_new, both_meta$Source,
                         control = how(nperm=999))
summary(indval_both, alpha = 1)

#plot of indicator species - supplementary figure
relAbdPlot <- list()
i <- 1
colours <- c( "#FDAE99","#882255","#DDCC77","#aa4499","#CC6677","#117733","#88ccee","#332288","#44aa99")
for(i in 1:9){
  coords_sp <- cbind(coords_both, both_df_new[,i])
  colnames(coords_sp) <- c('MDS1','MDS2','Source','rel_abd')
  coords_sp_plot <- coords_sp[which(coords_sp$rel_abd > 0),]
  cols <- colours[i]
  #plot nmds coloured by gradient
  p_sp = ggplot(coords_sp_plot, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_point(shape = 21, stroke = 1.0, colour = cols, aes(size = rel_abd)) +
    geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 1.0, colour = "#000000", alpha = 0.5) +
    geom_point(data = coords_nmds_model, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 1.0, colour = "#888888", alpha = 0.25) +
    coord_cartesian(xlim = c(-0.3,0.25), ylim = c(-0.2, 0.25)) +
    #scale_color_manual(values = colour_names) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')
  p_sp
  relAbdPlot[[length(relAbdPlot)+1]] <- p_sp
}

relAbdAll <- plot_grid(
  plotlist = relAbdPlot,
  ncol = 3,
  byrow = TRUE
)
relAbdAll

ggsave("~/Documents/for_paper/figures_dump/4_nmds_relsp.pdf", width=6.4, height=6.4, units = "in")

#-------- additional analyses and supplementary figures --------

#model with shuffled aijs
library(ggplot2)
library(cowplot)
library(grid)
library(dplyr)
library(vegan)
library(indicspecies)
library(tidyr)
#open the file 
new_data_main <- as.data.frame(read.csv("~/Documents/cluster/endpoints_aij_shuffled.csv", header = T, row.names = 1))
new_data_plot <- as.data.frame(read.csv("~/Documents/cluster/endpoints_aij_shuffled.csv", header = T))
new_data_main <- round(new_data_main)
new_data_plot[,c(2:10)] <- round(new_data_plot[,c(2:10)]) #to get rid of the one spurious negative value that somehow snuck into the data

new_data_main[new_data_main < 100000] = 0 #remove species below detection limit

colour_names <- c("F.sanfranciscensis" = "#FDAE99",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99")

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum",
                  "C.paralimentarius","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

#processing the data so that it is relative abundance instead of absolute abundance
new_data <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(new_data) <- colnames(new_data_main)
for(i in 1:100){
  temp_data <- new_data_main[i,]
  temp_data <- temp_data/sum(temp_data)
  new_data <- rbind(new_data,temp_data)
}

#community composition
#distance matrix - bray-curtis dissimilarity
distance_model <- vegdist(new_data, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_model <- hclust(distance_model, method = "ward.D2")
plot(hc_model)

dd_model <- as.dendrogram(hc_model)
pdf("~/Documents/for_paper/figures_dump/sup2_model_dend.pdf", width=13.6, height=3.20) # The height of the plot in inches
plot(dd_model) #1242X521
dev.off()

#make a data frame with the cluster order
model_clustorder<-hc_model$labels[c(hc_model$order)] 
model_clustorder<-as.data.frame(model_clustorder)
model_clustorder<-cbind(model_clustorder,rownames(model_clustorder))
colnames(model_clustorder) <- c("Population","order")
model_clusters <- model_clustorder

# change data from wide format to long format
colIDs <- colnames(new_data)
new_data_long <- new_data_plot %>% pivot_longer(cols=all_of(colIDs),
                                                names_to='Species',
                                                values_to='Abundance')
new_data_long <- merge(new_data_long, model_clusters, by = 'Population')

new_data_long$order<-as.numeric(new_data_long$order)
new_data_long<-arrange(new_data_long, order)

p_ordered <- ggplot(data = new_data_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = factor(Species, levels = strain_names))) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_ordered #1242X261
ggsave("~/Documents/for_paper/figures_dump/sup2_aij0_clustered.pdf", width=12.42, height=2.60, units = "in")

#plotting alpha estimates to see if they are bounded at 0
library(ggplot2)
library(dplyr)
library(cowplot)
#reading data from the pairwise parameters sheet - new data
paired_data <-  read.table(file = "/Users/kasturilele/Documents/community/est_all_pairs_6-3.txt", sep = ",", header = TRUE)

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum","C.paralimentarius","S.cerevisiae", "W.anomalus","K.humilis","K.servazzii")
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")
colour_strainID <- c("163" = "#88CCEE",
                     "228" = "#332288",
                     "253" = "#117733",
                     "177" = "#44AA99",
                     "17B2" = "#FDAE99",
                     "0092a" = "#882255",
                     "232" = "#DDCC77",
                     "460" = "#CC6677",
                     "550" = "#AA4499") #colour names with strain IDs
# #order data by correct order
paired_data$Strain_1 <- factor(paired_data$Strain_1, levels=strain_order)
paired_data$Strain_2 <- factor(paired_data$Strain_2, levels=strain_order)

#plot desnity plots of parameters
#plotting r (since they were not bounded by anything)
plot1 <- ggplot(data = paired_data) +
  theme_bw() +
  geom_density(aes(r1,colour = Strain_1)) +
  geom_density(aes(r2,colour = Strain_2)) +
  facet_grid(Strain_2~Strain_1) +
  scale_color_manual(values = colour_strainID) +
  labs(x = "growth rate", y = "count") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
plot1  

#plotting aii 
plot2 <- ggplot(data = paired_data) +
  theme_bw() +
  #geom_density(aes(log10(-a11),colour = Strain_1)) +
  #geom_density(aes(log10(-a22),colour = Strain_2)) +
  geom_density(aes(a11,colour = Strain_1)) +
  geom_density(aes(a22,colour = Strain_2)) +
  facet_grid(Strain_2~Strain_1) +
  scale_color_manual(values = colour_strainID) +
  labs(x = "aii", y = "count") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
plot2 

#plotting aij 
plot3 <- ggplot(data = paired_data) +
  theme_bw() +
  geom_density(aes(a12,colour = Strain_1)) +
  geom_density(aes(a21,colour = Strain_2)) +
  #geom_density(aes(log10(-a12),colour = Strain_1)) +
  #geom_density(aes(log10(-a21),colour = Strain_2)) +
  facet_grid(Strain_2~Strain_1, scales = "free") +
  scale_color_manual(values = colour_strainID) +
  labs(x = "aii", y = "count") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")
plot3 

#having each individual plot be on its own scale (do this simultaneously for all 3 plots)
count <- 1
plot_list_r <- list()
plot_list_aii <- list()
plot_list_aij <- list()
for(i in 1:8){
  S1 = strain_order[i]
  for(j in 2:9){
    if(j <= i){
      plot_list_r[[length(plot_list_r)+1]] <- pNew1
      plot_list_r[length(plot_list_r)] <- list(NULL)
      plot_list_aii[[length(plot_list_aii)+1]] <- pNew2
      plot_list_aii[length(plot_list_aii)] <- list(NULL)
      plot_list_aij[[length(plot_list_aij)+1]] <- pNew3
      plot_list_aij[length(plot_list_aij)] <- list(NULL)
      print(paste(length(plot_list_r),length(plot_list_aii),length(plot_list_aij)))
    }
    else{
      S2 = strain_order[j]
      cat(S1, ", ", S2, "\n")
      print(count)
      flip <- FALSE
      temp_pars <- subset(paired_data, Strain_1 == S1 & Strain_2 == S2)
      if(nrow(temp_pars) < 1){
        temp_pars <- subset(paired_data, Strain_1 == S2 & Strain_2 == S1)
        flip <- TRUE #S1 and S2 are flipped, this affects calculation of f2/f1
      }
      # pNew1 <- ggplot(data = temp_pars) +
      #   theme_bw() +
      #   geom_density(aes(r1,colour = Strain_1)) +
      #   geom_density(aes(r2,colour = Strain_2)) +
      #   scale_color_manual(values = colour_strainID) +
      #   labs(x = "r", y = "count") + 
      #   theme(axis.line = element_blank(),
      #         axis.text = element_text(size = 6),
      #         axis.ticks = element_line(linewidth = 0.1),
      #         legend.position="none")
      # plot_list_r[[length(plot_list_r)+1]] <- pNew1
      # pNew2 <- ggplot(data = temp_pars) +
      #   theme_bw() +
      #   geom_density(aes(a11,colour = Strain_1)) +
      #   geom_density(aes(a22,colour = Strain_2)) +
      #   scale_color_manual(values = colour_strainID) +
      #   labs(x = "aii", y = "count") + 
      #   theme(axis.line = element_blank(),
      #         axis.text = element_text(size = 6),
      #         axis.ticks = element_line(linewidth = 0.1),
      #         legend.position="none")
      # plot_list_aii[[length(plot_list_aii)+1]] <- pNew2
      pNew3 <- ggplot(data = temp_pars) +
        theme_bw() +
        geom_density(aes(a12,colour = Strain_1)) +
        geom_density(aes(a21,colour = Strain_2)) +
        scale_color_manual(values = colour_strainID) +
        labs(x = "aij", y = "count") + 
        theme(axis.line = element_blank(),
              axis.text = element_text(size = 6),
              axis.ticks = element_line(linewidth = 0.1),
              legend.position="none")
      plot_list_aij[[length(plot_list_aij)+1]] <- pNew3
      count = count + 1
    }
  }
}

#plot all the bar plots in a grid
p_all_r <- plot_grid(
  plotlist = plot_list_r,
  ncol = 8,
  byrow = TRUE
)
p_all_r #1048X1048
ggsave("~/Documents/for_paper/revisions/figdump/sup_all_r.pdf", width=12.00, height=10.48, units = "in")
#plot all the bar plots in a grid
p_all_aii <- plot_grid(
  plotlist = plot_list_aii,
  ncol = 8,
  byrow = TRUE
)
p_all_aii #1048X1048
ggsave("~/Documents/for_paper/revisions/figdump/sup_all_aii.pdf", width=14.00, height=10.48, units = "in")
#plot all the bar plots in a grid
p_all_aij <- plot_grid(
  plotlist = plot_list_aij,
  ncol = 8,
  byrow = TRUE
)
p_all_aij #1048X1048
ggsave("~/Documents/for_paper/revisions/figdump/sup_all_aij.pdf", width=14.00, height=10.48, units = "in")


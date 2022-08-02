#perform dbscan clustering & project clustering on tsne

# STEP 1: Importing Necessary Libraries

rm(list=ls())
library(Rtsne)
library(ggplot2)
require("ggrepel")
# For Data Manipulation
library(tidyverse) 
# For Clustering algorithm
library(cluster)
library(fpc)
library(mclust)
# for cluster visualisation
library(factoextra)

featurename <- ""

#STEP 2: Loading the Dataset
path <- "\\dir\\"
myData<-read.csv(paste0(path,"Top20.csv")  ,header=T,sep=",",quote=NULL, comment="")

#Renaming in full form
myData$Lithotype[myData$Lithotype=="Srp Per"]<- "Ultrabasic"
myData$Lithotype[myData$Lithotype=="Chr"]<- "Chromitite"
myData$Lithotype[myData$Lithotype=="Gbr"]<- "Basic"
myData$Lithotype[myData$Lithotype=="Rdg Gbr"]<- "Rodingtised Basic"
myData$Lithotype[myData$Lithotype=="Bst"]<- "Volcanic Basic"

str(myData)
colnames(myData)

numericData<-myData[,c(2:22)]
colnames(numericData)

numericData<- numericData[,c(1,9,16,5, 18, 21,
                             8,2,3,4, 14,13, 
                             15,20, 6, 19, 7,
                             17,12,11,10 )]
colnames(numericData)

# Select feature set by uncommenting the required three lines of code

# CH4,V & Fe
# numericData<- numericData[,c(1, 15, 16)]
# colnames(numericData)
# featurename <- "CH4VFe"

# V & Fe (noCh4)
# numericData<- numericData[,c(15, 16)]
# colnames(numericData)
# featurename <- "VFe_noCH4"

# #Top2_noCH4
# numericData<- numericData[,c(2:3)]
# colnames(numericData)
# featurename <- "Top2_noCH4"

# #top2
# numericData<- numericData[,c(1:3)]
# colnames(numericData)
# featurename <- "Top2"

# #Top3_noCH4
# numericData<- numericData[,c(2:4)]
# colnames(numericData)
# featurename <- "Top3_noCH4"

# #top3
# numericData<- numericData[,c(1:4)]
# colnames(numericData)
# featurename <- "Top3"


# # #Top4_noCH4
# numericData<- numericData[,c(2:5)]
# colnames(numericData)
# featurename <- "Top4_noCH4"

# #top4
# numericData<- numericData[,c(1:5)]
# colnames(numericData)
# featurename <- "Top4"

# #Top7_noCH4
# numericData<- numericData[,c(2:8)]
# colnames(numericData)
# featurename <- "Top7_noCH7"

#top7
# numericData<- numericData[,c(1:8)]
# colnames(numericData)
# featurename <- "Top7"
#

# #Top12_noCH4
# numericData<- numericData[,c(2:13)]
# colnames(numericData)
# featurename <- "Top12_noCH4"

#
# #top12
# numericData<- numericData[,c(1:13)]
# colnames(numericData)
# featurename <- "Top12"


# #top14
# numericData<- numericData[,c(1:15)]
# colnames(numericData)
# featurename <- "Top14"  #seed 111

# #Top14_noCH4
# numericData<- numericData[,c(2:15)]
# colnames(numericData)
# featurename <- "Top14_noCH4"  #seed 111

# 
# #top16
# numericData<- numericData[,c(1:17)]
# colnames(numericData)
# featurename <- "Top16"

# #top16
# numericData<- numericData[,c(2:17)]
# colnames(numericData)
# featurename <- "Top16_noCH4"

##CLUSTERING#############################################
#https://www.projectpro.io/recipes/do-dbscan-clustering-r
# getting the required information about the dataset
# glimpse(numericData)

#STEP 3: Data Preprocessing (Scaling)
# scaling the dataset
# numericData = scale(numericData)
# numericData %>% head()


# STEP 4: Obtaining Optimal value of eps
# We use the kNNdistplot(data, k=) function to carry this out task. It calculates the radius of the clusters.
# The method proposed here consists of computing the k-nearest neighbor distances in a matrix of points.
# The idea is to calculate, the average of the distances of every point to its k nearest neighbors. 
# The value of k will be specified by the user and corresponds to MinPts.
# Next, these k-distances are plotted in an ascending order. The aim is to determine the “knee”, which corresponds to the optimal eps parameter.
# A knee corresponds to a threshold where a sharp change occurs along the k-distance curve.

dir.create(paste0("dir\\",featurename))
set.seed(111)
BIC <- mclustBIC(numericData)
plot(BIC)

sink(file=paste0("dir\\",featurename,"\\",featurename,".txt"))

summary(BIC)
db <- Mclust(numericData, x = BIC, G=5)
# plot(db, what = "classification")
summary(db, parameters = TRUE)

table(db$classification,myData$Lithotype)

db
print("==========db$modelName==========")
db$modelName

print("==========db$BIC==========")
db$BIC

print("==========db$classification==========")
db$classification

print("==========db$uncertainty==========")
db$uncertainty

library(vcd)
print("=====CH4 Levels=====")
table(db$classification, myData$CH4label)
print(assocstats(xtabs(~as.factor(myData$CH4label)+db$classification)))
summary(xtabs(~as.factor(myData$CH4label)+db$classification))
summary(assocstats(xtabs(~as.numeric(as.factor(myData$CH4label))+db$classification)))
# chisq.test(xtabs(~as.numeric(as.factor(myData$CH4label))+db$classification), simulate.p.value = TRUE)
#It gave the warning because many of the expected values will be very small and therefore the approximations of p may not be right.
# "Assumptions" section of Pearson's chi-squared test article. In a nutshell, when expected (not actual) counts in any of the cells in your table 
# are fewer than 5 then one of the assumptions is broken, this is true for 1 DF situation; not >1 DF

print("====GoodmanKruskalTau====")
# https://cran.r-project.org/web/packages/GoodmanKruskal/vignettes/GoodmanKruskal.html
library(GoodmanKruskal)
GKtau(myData$CH4label,db$classification)

library(greybox)
assoc(myData$CH4label,db$classification)

print("=====Lithotype=====")
table(db$classification, myData$Lithotype)
print(assocstats(xtabs(~as.numeric(as.factor(myData$Lithotype))+db$classification, data=myData)))
summary(xtabs(~as.numeric(as.factor(myData$Lithotype))+db$classification, data=myData))
summary(assocstats(xtabs(~as.numeric(as.factor(myData$Lithotype))+db$classification, data=myData)))
# chisq.test(xtabs(~as.numeric(as.factor(myData$Lithotype))+db$classification, data=myData), simulate.p.value = TRUE)

clusnum<-5
class <- myData$Lithotype 

print("###CLUSTER#######################################")
for(i in 1:clusnum){
  print(paste0("::::cluster=",i,"::::"))
  print(myData$Sample[db$classification==i])
}
print("#####LITHO####################################################")
for(i in 1:clusnum){
  print(paste0("::::cluster=",i,"::::"))
  print(class[db$classification==i])
}


sink()

d <- data.frame(rbind(apply(numericData[db$classification==1,],2,mean), apply(numericData[db$classification==2,],2,mean)))
d<- rbind(d,apply(numericData[db$classification==3,],2,mean))
d<- rbind(d, apply(numericData[db$classification==4,],2,mean))
d<- rbind(d, apply(numericData[db$classification==5,],2,mean))

write.table(d,paste0("dir\\",featurename,"\\k",clusnum,"centres_",featurename,".csv") , sep="\t")

###DRAWING##########################################
type<-factor(db$classification)
clusnum<-length(unique(db$classification))
CH4label <- factor(myData$CH4label, levels = c("Very High", "High","Moderate", "Low", "Negligible"))
Litho <- factor(myData$Lithotype, levels = c("Chromitite", "Ultrabasic", "Basic", "Rodingtised Basic", "Volcanic Basic"))
# col_scheme <- c("blue", "red", "green", "purple", "magenta","cyan", "darkgoldenrod", "black", "orange", "darkgreen", "darkblue", 
#                 "darkred", "coral", "deeppink", "darkorange", "yellowgreen")

#safe_colorblind_palette 
# col_scheme<- c("#888888", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                "#44AA99", "#999933", "#882255", "#661100" )
##FULL COLOUR SCHEME ################################################
#PiYG =     rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef", 
#"#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221", "#276419"))
##################################################

col_scheme<- c("#8e0152", "#276419", "#de77ae", "#7fbc41", "#f1b6da", "#b8e186")
# scores_plot<-NULL
# newdf <- data.frame(myData$Ir, myData$Ru)
# colnames(newdf)<- c("Ir","Ru")
# str(newdf)
# scores_plot <- ggplot(newdf,  aes(x=Ir, y=Ru)) + theme_bw() + theme(panel.grid.major = element_blank(),
#                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# scores_plot<- scores_plot + ggtitle(paste0("Ru vs Ir",featurename)) + xlab("Ir") + ylab("Ru")
# scores_plot <- scores_plot + geom_point(aes(shape= Litho,fill=CH4label),size=4, color="white")
# scores_plot <- scores_plot +  theme(aspect.ratio=1)
# scores_plot <- scores_plot + scale_colour_manual(name="Cluster", values=col_scheme) # for clusters
# scores_plot <- scores_plot + scale_fill_manual(name="CH4 Level", values=col_scheme[1:5])
# scores_plot <- scores_plot + scale_shape_manual(name="Lithotype", values=c(21,22,23,24,25))
# 
# scores_plot <- scores_plot + geom_text(check_overlap = FALSE, hjust = "inward",vjust="inward", nudge_x = 0.01,
#                                        aes(label = myData$Sample, colour = type), size=3) #for textcolour of clusters
# 
# scores_plot <- scores_plot + guides(colour =guide_legend(override.aes = list(colour=col_scheme[1:clusnum], size=4)),
#                                     fill =guide_legend(override.aes = list(colour=col_scheme[1:5], size=4)),
#                                     shape =guide_legend(override.aes = list(shape=c(21,22,23,24,25),colour="white",
#                                                                             fill=c(1,1,1,1,1),size=4)))
# 
# # scores_plot <- scores_plot + geom_mark_ellipse(expand=0, color="grey", aes(label=tsne, group=tsne),
# #                                                label.colour=4, label.buffer = unit(10,"mm"), 
# #                                                con.type="straight", con.colour = "grey") #, aes(fill=periods)
# png(filename=paste0("dir\\top",dim(numericData)[2]-1, "_IrRu.png"), width = 650, height = 650)
# #when in a loop, put print
# print(scores_plot)
# dev.off()

#set.seed(155)
#seed<-140
#set.seed(seed)
# plex<-19
# for(i in 1:5){
for(plex in 1:19){
  
  print(plex)
  maxiter<-6000
  tsne_results <- tryCatch(Rtsne(numericData, perplexity=plex, 
                                 check_duplicates = FALSE,
                                 max.iter=maxiter),
                           error=function(cond) {
                             message("Error with tsne result")
                             return(-1)
                           } )# You can change the value of perplexity and see how the plot changes
  
  if (length(tsne_results) <2) next
  # print(paste0("success=", plex))
  # print(tsne_results)
  
  #####################################################################################################
  tsnedata<- data.frame(tsne_results$Y)
  scores_plot<-NULL
  scores_plot <- ggplot(tsnedata,  aes(x=X1, y=X2)) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  scores_plot<- scores_plot + ggtitle(paste0("t-SNE",featurename,"_plex=",plex,
                                             "_maxIter=", maxiter)) + xlab("X1") + ylab("X2")
  scores_plot <- scores_plot + geom_point(aes(shape= Litho,fill=CH4label),size=4, color="white")
  scores_plot <- scores_plot +  theme(aspect.ratio=1)
  scores_plot <- scores_plot + scale_colour_manual(name="Cluster", values=col_scheme) # for clusters
  scores_plot <- scores_plot + scale_fill_manual(name="CH4 Level", values=col_scheme[1:5])
  scores_plot <- scores_plot + scale_shape_manual(name="Lithotype", values=c(21,22,23,24,25))
  
  scores_plot <- scores_plot + geom_text(check_overlap = FALSE, hjust = "inward",vjust="inward", nudge_x = 0.01,
                                         aes(label = myData$Sample, colour = type), size=3) #for textcolour of clusters
  
  scores_plot <- scores_plot + guides(colour =guide_legend(override.aes = list(colour=col_scheme[1:clusnum], size=4)),
                                      fill =guide_legend(override.aes = list(colour=col_scheme[1:5], size=4)),
                                      shape =guide_legend(override.aes = list(shape=c(21,22,23,24,25),colour="white",
                                                                              fill=c(1,1,1,1,1),size=4)))
  
  # scores_plot <- scores_plot + geom_mark_ellipse(expand=0, color="grey", aes(label=tsne, group=tsne),
  #                                                label.colour=4, label.buffer = unit(10,"mm"), 
  #                                                con.type="straight", con.colour = "grey") #, aes(fill=periods)
  png(filename=paste0("dir\\",featurename,"\\",featurename,"_plex=",plex, ".png"), width = 650, height = 650)
  #when in a loop, put print
  print(scores_plot)
  dev.off()
  
  #####################################################################################################

  # }
}
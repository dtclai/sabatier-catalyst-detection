#t-sne visualisation
rm(list=ls())
library(Rtsne)
library(ggplot2)
require("ggrepel")

featurename <- ""


path <- "~\\dir\\"
myData<-read.csv(paste0(path,"Top20.csv")  ,header=T,sep=",",quote=NULL, comment="")

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
# featurename <- "Top7_noCH4"

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

CH4label <- factor(myData$CH4label, levels = c("Very High", "High","Moderate", "Low", "Negligible"))
Litho <- factor(myData$Lithotype, levels = c("Chromitite", "Ultrabasic", "Basic", "Rodingtised Basic", "Volcanic Basic"))
tsne <- factor(myData$tsne)
# col_scheme <- c("blue", "red", "green", "purple", "magenta","cyan", "darkgoldenrod", "black", "orange", "darkgreen", "darkblue", 
#                 "darkred", "coral", "deeppink", "darkorange", "yellowgreen")

#safe_colorblind_palette 
# col_scheme<- c("#888888", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                "#44AA99", "#999933", "#882255", "#661100" )
# col_scheme<- c("#8e0152", "#276419", "#c51b7d", "#7fbc41", "#de77ae", "#b8e186")
##FULL COLOUR SCHEME ################################################
#PiYG =     rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef", 
#"#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221", "#276419"))
##################################################

col_scheme<- c("#8e0152", "#276419", "#de77ae", "#7fbc41", "#f1b6da", "#b8e186")

tsne <- factor(myData$tsne)
dir.create(paste0("~\\dir\\",featurename))


#set.seed(155)
#seed<-140
#set.seed(seed)
# plex<-19
# for(i in 1:5){
for(plex in 2:19){
 
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
  print(paste0("success=", plex))
  print(tsne_results)
  
  #####################################################################################################
  tsnedata<- data.frame(tsne_results$Y)
  scores_plot<-NULL
  scores_plot <- ggplot(tsnedata,  aes(x=X1, y=X2)) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  scores_plot<- scores_plot + ggtitle(paste0("t-SNE",featurename,"_plex=",plex,
                      "_maxIter=", maxiter)) + xlab("X1") + ylab("X2")
  scores_plot <- scores_plot + geom_point(aes(shape= Litho,fill=CH4label),size=4, color="white")
  scores_plot <- scores_plot +  theme(aspect.ratio=1)
  # scores_plot <- scores_plot + scale_colour_manual(name="Cluster", values=col_scheme) #no clusters
  scores_plot <- scores_plot + scale_fill_manual(name="CH4 Level", values=col_scheme[1:5])
  scores_plot <- scores_plot + scale_shape_manual(name="Lithotype", values=c(21,22,23,24,25))
  scores_plot <- scores_plot + geom_text(check_overlap = FALSE, hjust = "inward",vjust="inward", nudge_x = 0.01,
                                         aes(label = myData$Sample), size=3)
  # scores_plot <- scores_plot + geom_text(check_overlap = FALSE, hjust = "inward",vjust="inward", nudge_x = 0.01,
  #                                        aes(label = myData$Sample, colour = type), size=3) #no colour for clusters
  
  scores_plot <- scores_plot + guides(colour =guide_legend(override.aes = list(colour=col_scheme[1:5], size=4)),
                                      fill =guide_legend(override.aes = list(colour=col_scheme[1:5], size=4)),
                                      shape =guide_legend(override.aes = list(shape=c(21,22,23,24,25),colour="white",
                                                                              fill=c(1,1,1,1,1),size=4)))
  
  # scores_plot <- scores_plot + geom_mark_ellipse(expand=0, color="grey", aes(label=tsne, group=tsne),
  #                                                label.colour=4, label.buffer = unit(10,"mm"), 
  #                                                con.type="straight", con.colour = "grey") #, aes(fill=periods)
  
  png(filename=paste0("~\\dir\\",featurename,"\\",featurename,"_plex=",plex, ".png"), width = 650, height = 650)
  #when in a loop, put print
  print(scores_plot)
  dev.off()
# }
 }
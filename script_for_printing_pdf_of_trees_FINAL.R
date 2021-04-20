# change to directory
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
#setwd("C://Users/Edeline/Downloads/SANGER/04_21_2020/individual")
#setwd("C://Users/Edeline/Downloads/SANGER/04_27_2020")
setwd("E://Solanum_Phylogeny_Project/02_supermatrix/02_raxml/06_23_2020/00_RESULTS")

#ADD CLADE2 names to tree.
# load th file that contains the clade info.
clade<-read.csv("solanum_clades_acc_mod.csv",header=TRUE,sep=";")
colnames(clade)
dim(clade)
colnames(clade)[1]<-"FULLNAME"

species<-word(clade$FULLNAME,2,sep=" ")
species
clade$newname<-paste("S._",species,"_",clade$CLADE2,sep="")
clade$oldname<-paste("S._",species,sep="")

head(clade)
clade$newname[1]<-"J._procumbens_Jaltomata"
clade$newname[2]<-"J._sinuosa_Jaltomata"
clade$newname[3]<-"J._dentata_Jaltomata"
clade$newname[4]<-"J._bicolor_Jaltomata"

clade$oldname[1]<-"J._procumbens"
clade$oldname[2]<-"J._sinuosa"
clade$oldname[3]<-"J._dentata"
clade$oldname[4]<-"J._bicolor"



n<-c("Solanum paposanum", "Potato", "Regmandra", "S._paposanum_Regmandra", "S._paposum")
rbind(clade,n)->clade

n<-c("Solanum dimorphandrum", "Thelopodium", "Thelopodium", "S._dimorphandrum_Thelopodium", "S._dimorphamdrum")
rbind(clade,n)->clade

n<-c("Solanum mariae", "Potato", "Basarthrum", "S._mariae_Basarthrum", "S._mariaea")
rbind(clade,n)->clade

n<-c("Solanum nigrescens", "Morelloid", "Morelloid_Sect_Solanum", "S._nigrescens_Morelloid_Sect_Solanum", "S._nigriscens")
rbind(clade,n)->clade

n<-c("Solanum salamancae", "Morelloid", "Morelloid_Sect_Solanum", "S._salamancae_Morelloid_Sect_Solanum", "S._salamanceae")
rbind(clade,n)->clade

n<-c("Solanum xblanco-galdosii", "Potato", "Petota", "S._blanco???galdosii_Petota", "S._blanco???galdosii")
rbind(clade,n)->clade

n<-c("Solanum plowmanii", "Geminata", "Geminata", "S._plowmanii_Geminata", "S._plowmannii")
rbind(clade,n)->clade

n<-c("Solanum polhillii", "Leptostemonum", "Old_World-Africa", "S._polhillii_Old_World-Africa", "S._polhilii")
rbind(clade,n)->clade

n<-c("Solanum velardei", "Potato", "Petota", "S._velardei_Petota", "S._verladei")
rbind(clade,n)->clade

n<-c("Solanum glutinosum", "Leptostemonum", "Torva", "S._glutinosum_Torva", "S._glutinosa")
rbind(clade,n)->clade

n<-c("Solanum glutinosum", "Leptostemonum", "Torva", "S._glutinosum_Torva", "S._glutinosa")
rbind(clade,n)->clade

n<-c("Solanum knoblochii", "Leptostemonum", "Androceras", "S._knoblochii_Torva", "S._knoblichi")
rbind(clade,n)->clade



recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}


##########

f.tr<-dir(pattern=".T15")
f.tr

#pdf(width=8.5, height=40,file="Solanum_results_individual_DATA_June_2020_V2.pdf")

#list<-c()

for (i in 1:2)
{
i<-3
  print(f.tr[i])
  tr<-treeio::read.tree(file.choose())
  #This is for the final tree.
#  tr<-drop.tip(tr,c("S._nigriscens","S._viride","S._probolospermum","S._urubambaense"))
  
  #search and replace the tip names of the trees. 
  
  newnames<-recoderFunc(tr$tip.label, clade$oldname, clade$newname)
  tr$tip.label<-newnames
  tr2<-tr
  tr2<-root(tr,"J._bicolor_Jaltomata")
  
  
  nodelabel<-data.frame(nodelabel=tr2$node.label, node=tr2$Nnode)
  #colnames(nodelabel)<-c("parent", "node", "nodelabel")
  
  q <- ggtree(tr2)
  d <- q$data
  d <- d[!d$isTip,]
  d$label <- as.numeric(d$label)
  d <- d[d$label > 50,]
  d
  
  z<-ggtree(tr2, ladderize=T)+
    
    #Add a title
    ggtitle(f.tr[i])+
    # Adds tip labels
    geom_tiplab(color="black",size=1.5) +
    geom_text2(aes(subset = !isTip, label=label),size=1.5, hjust=1.25,vjust=1.5,fontface="italic")
  
  z + geom_text(data=d, aes(label=label))
  sanger.z<-z + geom_point2(data=d, aes(subset=label>=75),size=2, shape=21, fill="red") + geom_point2(data=d, aes(subset=label>=95),size=2, shape=21, fill="cyan")
  ggsave(sanger.z, file=paste(f.tr[i],"2.pdf",sep=""), width=8.5, height=35)
  
}

plot(list[1])
dev.off()


# list all result files in directory
f.tr<-dir(pattern="*.T15")
f.tr

pdf(width=8.5, height=40,file="Solanum_results_individual_DATA_25_May_2020.pdf")

for (i in 1:length(f.tr))
{
  print(f.tr[i])
  tr<-read.tree(f.tr[i])
  tr2<-root(tr,"J._bicolor")
  
  #search and replace the tip names of the trees. 
  
  newnames<-recoderFunc(tr2$tip.label, clade$oldname, clade$newname)
  tr2$tip.label<-newnames
  tr2
  
  plot(tr2, cex=0.4)
  nodelabels(tr2$node.label,cex=0.4,frame="n",adj=c(1,0))
  title(f.tr[i],cex=0.5)
}

dev.off()


##########################
f.tr<-dir(pattern="*.T15")
f.tr
# For each 

pdf(width=8.5, height=45,file="Solanum_results_ALL_DATA_25_May_2020.pdf")

for (i in 1:length(f.tr))
{
  print(f.tr[i])
  tr<-read.tree(f.tr[i])
  tr2<-root(tr,"J._bicolor")
  
  #search and replace the tip names of the trees. 
  
  newnames3<-recoderFunc(tr2$tip.label, clade$oldname, clade$newname)
  tr2$tip.label<-newnames
  tr2
  
  plot(tr2, cex=0.4)
  nodelabels(tr2$node.label,cex=0.4,frame="n",adj=c(1,0))
  title(f.tr[i],cex=0.5)
}

dev.off()

##############################


# list all result files in directory
f.tr<-dir(pattern="*.T15")
f.tr
# For each 

pdf(width=8.5, height=45,file="Solanum_results_ALL_DATA_JUNE_2020_V2.pdf")

for (i in 1:length(f.tr))
{
  
  print(f.tr[i])
  tr<-read.tree(f.tr[i])
  
  #search and replace the tip names of the trees. 
  
  newnames<-recoderFunc(tr$tip.label, clade$oldname, clade$newname)
  tr$tip.label<-newnames
  tr
  tr2<-root(tr,"J._bicolor_Jaltomata")
  
  #nodelabel<-data.frame(nodelabel=tr2$node.label, node=tr2$Nnode)
  
  #colnames(nodelabel)<-c("parent", "node", "nodelabel")
  
  q <- ggtree(tr2)
  d <- q$data
  d <- d[!d$isTip,]
  d$label <- as.numeric(d$label)
  d <- d[d$label > 50,]
  d
  
  z<-ggtree(tr2, ladderize=T)+
    
    #Add a title
    ggtitle(f.tr[i])+
    # Adds tip labels
    geom_tiplab(color="black",size=1.5) +
    geom_text2(aes(subset = !isTip, label=label),size=1.5, hjust=1.25,vjust=1.5,fontface="italic")
  
  z + geom_text(data=d, aes(label=label))
  sanger.z<-z + geom_point2(data=d, aes(subset=label>75),size=2, shape=21, fill="orange") + geom_point2(data=d, aes(subset=label>94),size=2, shape=21, fill="chartreuse")
  sanger.z  
}

dev.off()

###
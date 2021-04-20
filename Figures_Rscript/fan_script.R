
# CLEAN WORKSPACE ====================
rm(list=ls())

## DIRECTORIES ====================
main.dir <- "E://Solanum_Phylogeny_Project/02_supermatrix/02_raxml/06_23_2020/00_RESULTS"
fig.dir <- file.path(main.dir, "figs")
#data.dir <- file.path(main.dir, "data")

## LOAD PACKAGES ====================
# You will need to install ggtree which is crucial for the plotting in this script, run the following code to ensure that all dependent packages are also installed and updated
# install.packages(c("devtools", "rgdal", "rgeos", "ape", "plyr", "wesanderson", "cowplot", "scatterpie", "ggrepel", "ggplot2", "BiocManager", "viridis", "reshape2"), dependencies = TRUE)
# devtools::install_github('GuangchuangYu/ggtree', force = TRUE)

# spatial packages for maps
#library(rgdal); library(rgeos)

# tidy packages for data handling and summary statistics
library(plyr); library(reshape2)

# phylogenetic package
library(ape)

# general plotting packages
library(ggplot2);library(wesanderson); library(cowplot); library(ggrepel); library(scatterpie); library(ggtree); library(viridis)

# Load phylogeny
#     This is a Maximum Clade Credibility (MCC) phylogenetic tree as used by Onstein et al. (2017)
#     Nature Ecology & Evolution 1: 1903-1911.
#     and by Onstein et al. (2018), Proceedings of the Royal Society B: Biological Sciences, 285,
#     20180882.
#     Available from DRYAD: https://datadryad.org/resource/doi:10.5061/dryad.cm4nm
solPhylo <- read.tree(file.path(main.dir, "RAxML_allData_v9.T15"))

#update.names
#####################
#ADD CLADE2 names to tree.
# load th file that contains the clade info.
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

n<-c("Solanum xblanco-galdosii", "Potato", "Petota", "S._blanco-galdosii_Petota", "S._blanco-galdosii")
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

solPhylo<- read.tree(file.path(main.dir, "RAxML_allData_v9.T15"))

solPhylo$tip.label <- gsub(solPhylo$tip.label, pattern = "_", replacement = " ")
solPhylo$tip.label


##########

newnames<-recoderFunc(solPhylo$tip.label, clade$oldname, clade$newname)
solPhylo$tip.label<-newnames
solPhylo$tip.label


########################################################################
## Fan tree

#Parameter in tree

bs = 2; bs1=1; fs = 3.5; ofs = 30; cl = c("grey90", "grey20")


phyloCladePlot <- ggtree(solPhylo, layout = "circular", size=0.2)+

   #geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 1) +
  # geom_tiplab2(size = 1)
 #xlim(c(-20, 10)) +
  geom_point2(aes(subset=!isTip&label>74, label=node),size=0.8, shape=21, fill="orange") + geom_point2(aes(subset=!isTip&label>=95, label=node),size=0.8, shape=21, fill="chartreuse")+
 
  #Potato
  #  MRCA(solPhylo,"S._okadae_Petota", "S._phaseoloides_Herpystichum")[1381]
  geom_cladelabel(node = 1381, label = "Potato", hjust = 0, offset= 0.04,offset.text = 0.01, barsize = bs, fontsize = 5, color = "blue") +
  
  #Petota
  #  MRCA(solPhylo,"S._okadae_Petota", "S._chomatophilum_Petota")[1437]
  geom_cladelabel(node = 1437, label = "Petato", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +
  
  #Tomato
  #  MRCA(solPhylo,"S._galapagense_Tomato", "S._juglandifolium_Tomato")[1417]
  geom_cladelabel(node = 1417, label = "Tomato", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +

  #Etuberosum
  #  MRCA(solPhylo,"S._etuberosum_Etuberosum", "S._palustre_Etuberosum")[1413]
  geom_cladelabel(node = 1413, label = "Etuberosum", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +
  
  #Anarrhichomenum
  #  MRCA(solPhylo,"S._complectens_Anarrhichomenum", "S._ionidium_Anarrhichomenum")[1405]
  geom_cladelabel(node = 1405, label = "Anarrhichomenum", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +
  
  #Basarthrum
  #  MRCA(solPhylo,"S._suaveolens_Basarthrum", "S._cochoae_Basarthrum")[1489]
  geom_cladelabel(node = 1489, label = "Basarthrum", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +
  
  #Articulatum
  #Oxycocoides
  
  #Pteroidea - Herpystichum
  #  MRCA(solPhylo,"S._anceps_Pteroidea", "S._phaseoloides_Herpystichum")[1382]
  geom_cladelabel(node = 1382, label = "Pteroidea - Herpystichum", hjust = 0, offset.text = 0.01, barsize = bs1, fontsize = fs, color = "blue") +
  
  
  #Regmandra
  #  MRCA(solPhylo,"S._multifidum_Regmandra", "S._brachyantherum_Regmandra")[1376]
  geom_cladelabel(node = 1376, label = "Regmandra", hjust = 0, offset= 0.05, offset.text = 0.01, barsize = bs, fontsize = 5, color = "blue") +

  #M Clade  
  #  MRCA(solPhylo,"S._endoadenium_Dulcamaroid", "S._macrothyrsum_African non-spiny")[1267]
  geom_cladelabel(node = 1267, label = "M Clade", hjust = 0, offset= 0.05, barsize = bs, fontsize = 5, color = "purple") +
  
  #Dulcamaroid + Morelloid  
  #  MRCA(solPhylo,"S._endoadenium_Dulcamaroid", "S._cochabambense_Morelloid_Sect_Solanum")[1402]
  #geom_cladelabel(node = 1402, label = "DulMo", hjust = 0, offset= 0.1, offset.text = 0.05, barsize = bs, fontsize = fs, color = "purple") +

  #Dulcamaroid
  #  MRCA(solPhylo,"S._endoadenium_Dulcamaroid", "S._lyratum_Dulcamaroid")[1285]
  geom_cladelabel(node = 1285, label = "Dulcamaroid", hjust = 0, offset.text = 0.05, barsize = bs1, fontsize = fs, color = "purple") +
  
  #Morelloid
  #  MRCA(solPhylo,"S._tripartitum_Morelloid_Radicans", "S._cochabambense_Morelloid_Sect_Solanum")[1309]
  geom_cladelabel(node = 1309, label = "Morelloid", hjust = 0, offset.text = 0.05, barsize = bs1, fontsize = fs, color = "purple") +
  
  #Vans clade
  #  MRCA(solPhylo,"S._laciniatum_Archaesolanum", "S._macrothyrsum_African non-spiny")[1268]
  geom_cladelabel(node = 1268, label = "Vans clade", hjust = 0, offset.text = 0.05, barsize = bs1, fontsize = fs, color = "purple") +
  
  #Thelopodium  
  #  MRCA(solPhylo,"S._dimorphandrum_Thelopodium", "S._monarchostemon_Thelopodium")[747]
 #geom_cladelabel(node = 753, label = "Thelopodium", hjust = 0, ofs= 1, offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +

  #Clade 2
  #  MRCA(solPhylo,"S._anomalostemon_Mapiriense", "S._parvifolium_Old_World-Australia")[758]
  geom_cladelabel(node = 758, label = "Clade II", hjust = 0, offset=0.05, offset.text = 0.05, barsize = bs, fontsize = 5, color = cl) +
  
  #Leptostemonum S._polygamum_Leptostemonum_unplaced S._parvifolium_Old_World-Australia
  #  MRCA(solPhylo,"S._polygamum_Leptostemonum_unplaced", "S._parvifolium_Old_World-Australia")[810]
  geom_cladelabel(node = 810, label = "Leptostemonum", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +
    
    #OldWorld S._virginianum_Old_World-Tropical Asia S._parvifolium_Old_World-Australia
    #  MRCA(solPhylo,"S._virginianum_Old_World-Tropical Asia", "S._parvifolium_Old_World-Australia")[874]
    geom_cladelabel(node = 874, label = "Old World", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
    #Torva S._peikuoense_Torva S._donianum_Torva
    #  MRCA(solPhylo,"S._peikuoense_Torva", "S._donianum_Torva")[1075]
    geom_cladelabel(node = 1075, label = "Torva", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Elaegnifolium S._houstonii_Elaeagnifolium S._mortonii_Elaeagnifolium
  #  MRCA(solPhylo,"S._houstonii_Elaeagnifolium", "S._mortonii_Elaeagnifolium")[870]
  geom_cladelabel(node = 870, label = "Elaegnifolium", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Multispinum/Asterophorum grade
  
  #Micrantha
  #  MRCA(solPhylo,"S._lanceifolium_Micracantha", "S._volubile_Micracantha")[1109]
  geom_cladelabel(node = 1109, label = "Micracantha", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Bahamense
  #  MRCA(solPhylo,"S._bahamense_Bahamense", "S._ensifolium_Bahamense")[1119]
  geom_cladelabel(node = 1119, label = "Bahamense", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Crotonoides
    #  MRCA(solPhylo,"S._crotonoides_Crotonoides", "S._woodburyi_Crotonoides")[1117]
  geom_cladelabel(node = 1117, label = "Crotonoides", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
    
  #Caroliniense
  #  MRCA(solPhylo,"S._perplexum_Carolinense", "S._comptum_Carolinense")[1121]
  geom_cladelabel(node = 1121, label = "Caroliniense", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Campechiense
  
  #Androceras S._tenuipes_Androceras S._angustifolium_Androceras
  #  MRCA(solPhylo,"S._tenuipes_Androceras", "S._angustifolium_Androceras")[851]
  geom_cladelabel(node = 851, label = "Androceras", hjust = 0, offset=0.01,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Crinitum S._mitlense_Crinitum S._urticans_Crinitum
  #  MRCA(solPhylo,"S._mitlense_Crinitum", "S._urticans_Crinitum")[842]
  geom_cladelabel(node = 842, label = "Crinitum", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +
  
  #Sisymbriifolium
  #  MRCA(solPhylo,"S._sisymbriifolium_Sisymbriifolium", "S._hasslerianum_Sisymbriifolium")[840]
  geom_cladelabel(node = 840, label = "Sisymbriifolium", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Erythrotrichum S._robustum_Erythrotrichum S._rhytidoandrum_Erythrotrichum
  #  MRCA(solPhylo,"S._robustum_Erythrotrichum", "S._rhytidoandrum_Erythrotrichum")[1142]
  geom_cladelabel(node = 1142, label = "Erythrotrichum", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Thomasiifolium S._buddleifolium_Thomasiifolium S._rupincola_Thomasiifolium
  #  MRCA(solPhylo,"S._buddleifolium_Thomasiifolium", "S._rupincola_Thomasiifolium")[1139]
  geom_cladelabel(node = 1139, label = "Thomasiifolium", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Gardneri S._microphyllum_Gardneri S._stenandrum_Gardneri
  #  MRCA(solPhylo,"S._microphyllum_Gardneri", "S._stenandrum_Gardneri")[1131]
  geom_cladelabel(node = 1131, label = "Gardneri", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Acanthophora S._mammosum_Acanthophora S._tenuispinum_Acanthophora
  #  MRCA(solPhylo,"S._mammosum_Acanthophora", "S._tenuispinum_Acanthophora")[813]
  geom_cladelabel(node = 813, label = "Acanthophora", hjust = 0, offset=0.03,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Lasiocarpa S._stramoniifolium_Lasiocarpa S._parvifolium_Old_World-Australia
  #  MRCA(solPhylo,"S._stramoniifolium_Lasiocarpa", "S._repandum_Lasiocarpa")[825]
  geom_cladelabel(node = 825, label = "Lasiocarpa", hjust = 0,offset.text = 0.05, barsize = bs1, fontsize = fs, color = cl) +
  
  #Geminata acuminatum_Geminata
  # MRCA(solPhylo,"S._acuminatum_Geminata", "S._delitescens_Reductum")[1158]
  geom_cladelabel(node = 1158, label = "Geminata", hjust = 0, offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +
  
  #Brevantehrum
  #MRCA(solPhylo,"S._goodspeedii_Brevantherum s.s.", "S._trachytrichium_Brevantherum s.s.")[1227]
  geom_cladelabel(node = 1227, label = "Brevantherum", hjust = 0, offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +
  
  #Cyphomandra
  #MRCA(solPhylo,"S._graveolens_Wendlandii-Allophyllum", "S._proteanthum_Pachyphylla")[762]
  geom_cladelabel(node = 762, label = "Cyphomandra", hjust = 0, offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +
  
  #Nemorense
  #MRCA(solPhylo,"S._reptans_Nemorense", "S._nemorense_Nemorense")[1155]
  geom_cladelabel(node = 1155, label = "Nemorense", hjust = 0, offset.text = 0.1, barsize = bs, fontsize = fs, color = cl) +

  #Wendlandii-Allophylum
  #MRCA(solPhylo,"S._refractum_Wendlandii-Allophyllum", "S._allophyllum_Wendlandii-Allophyllum")[803]
  geom_cladelabel(node = 803, label = "Wendlandii-Allophyllum", hjust = 0, offset.text = 0.05, barsize = bs, fontsize = fs, color = cl) +

#  MRCA(solPhylo,"S._anomalostemon_Mapiriense", "S._parvifolium_Old_World-Australia")[754]
#  geom_cladelabel(node = 1262, label = "Regmandra", hjust = 0, ofs= 1, offset.text = 0.05, barsize = bs, fontsize = fs, color = "purple") +
  

  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"))+ 
  geom_treescale(color='red')

phyloCladePlot

setwd("E:/Solanum_Phylogeny_Project/03_plastome/")

name<-"fan_script_Raxml_right"
ggsave("fan_script_Raxml_rightv2.pdf", phyloCladePlot, width = 6, height = 5.6)
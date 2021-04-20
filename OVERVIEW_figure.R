#Open libraries
library(ape)
library(stringr)
library(stringi)
library(treeio)
library(ggtree)
library(ggplot2)

#Navigate to Path
setwd("E://Solanum_Phylogeny_Project/02_supermatrix/02_raxml/06_23_2020/00_RESULTS")

#Modify names
clade<-read.csv("E://Solanum_Phylogeny_Project/02_supermatrix/02_raxml/06_23_2020/00_RESULTS/solanum_clades_acc_mod.csv",header=TRUE,sep=";")
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
#Read raxml file
tr3<-ape::read.tree("RAxML_allData_v9.T15")
tr5<-ape::read.tree("RAxML_bestTree.result")
#search and replace the tip names of the trees. 

newnames<-recoderFunc(tr3$tip.label, clade$oldname, clade$newname)
tr3$tip.label<-newnames
tr4<-root(tr3,"J._bicolor_Jaltomata")


nodelabel<-data.frame(nodelabel=tr4$node.label, node=tr4$Nnode)
#colnames(nodelabel)<-c("parent", "node", "nodelabel")

q <- ggtree(tr4)
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 50,]
d

z<-ggtree(tr4, ladderize=T)+
  
  #Add a title
  ggtitle("all_DATA_RAXML")+
  # Adds tip labels
  geom_tiplab(color="black",size=1.5) +
  geom_text2(aes(subset = !isTip, label=label),size=1.5, hjust=1.25,vjust=1.5,fontface="italic")

z + geom_text(data=d, aes(label=label))
sanger.z<-z + geom_point2(data=d, aes(subset=label>=70),size=2, shape=21, fill="yellow") + geom_point2(data=d, aes(subset=label>=80),size=2, shape=21, fill="orange") + geom_point2(data=d, aes(subset=label>=95),size=2, shape=21, fill="chartreuse")
ggsave(sanger.z, file="All_data_test.pdf", width=8.5, height=40)



#Write as a Fan.

f<-ggtree(tr4,  branch.length='none', layout='circular')+
  
  #Add a title
  #ggtitle("all_DATA_RAXML")+
  
  # Adds tip labels
  #geom_tiplab(color="black",size=1.5) +
  geom_text2(aes(subset = !isTip & label>=70, label=label),size=1.5, hjust=1.25,vjust=1.5,fontface="italic")

#z + geom_text(data=d, aes(label=label))
f + geom_point2(aes(subset=!isTip & label>=70),size=2, shape=21, fill="yellow") + geom_point2(aes(subset=!isTip & label>=80),size=2, shape=21, fill="orange") + geom_point2(aes(subset=!isTip & label>=95),size=2, shape=21, fill="chartreuse")
fan_phylo<-f + geom_point2(aes(subset=!isTip & label>=80),size=2, shape=21, fill="orange") + geom_point2(aes(subset=!isTip & label>=95),size=2, shape=21, fill="chartreuse")

ggsave(fan_phylo, file="All_data_circle_test.pdf", width=8.5, height=40)

#Labelling Leptostemomum
#860
fan_phylo + geom_cladelabel(node=860, label="Leptostemonum", angle=25)

fan_phylo + geom_cladelabel(node=860, label="Leptostemonum", angle=25)

###################
#Trying to get thicker branchnodes....
#https://groups.google.com/forum/#!topic/bioc-ggtree/1TRaS5To5vI
x=treeio::read.newick('RAxML_All_bipartitions.T15.RESULT',node.label='support')
edge = x@phylo$edge
edge = as.data.frame(edge)
edge = edge[edge[,2] <= Ntip(x), ]

df = x@data
merge(df, edge, by.x="node", by.y="V1", all.x=F, all.y=T) -> df2
df3= df2[, c(3, 2)]
colnames(df3) = colnames(df)
x@data = rbind(df3, df)

ggtree(x, aes(size=I(support/100)))+ geom_point2(aes(subset=support==100),size=2, shape=21, fill="chartreuse")+
  geom_text2(aes(subset = !isTip & support>=70, label=support),size=1.5, hjust=1.25,vjust=1.5,fontface="italic")


###############

library(ggtree)
p = ggtree(tr4, ladderize=T) + geom_tiplab(size=1.5) + ggtitle("ALL_data_Raxml_NODE")
edge=data.frame(tr4$edge, edge_num=1:nrow(tr4$edge))
colnames(edge)=c("parent", "node", "edge_num")
p %<+% edge + geom_label(aes(label=node), size=1.5)
p %<+% edge + geom_text2(aes(label=node),size=1.5, hjust=1.25,vjust=1.5)

p %<+% edge + geom_text2(aes(label=node),size=1.5, hjust=1.25,vjust=1.5)




######
#Outgroup
MRCA(tr4,"J._bicolor_Jaltomata","J._procumbens_Jaltomata")#751

#Thelopodium
MRCA(tr4,"S._thelopodium_Thelopodium","S._dimorphandrum_Thelopodium")#753

#Label Clade 1
MRCA(tr4,"S._valdiviense_Valdiviense", "S._sitiens_Tomato")#1266

#Label Clade 2
MRCA(tr4,"S._graveolens_Wendlandii-Allophyllum", "S._anomalostemon_Mapiriense")#758


####
#Label Major clades
#M-Clade
MRCA(tr4,"S._valdiviense_Valdiviense", "S._americanum_Morelloid_Sect_Solanum")

  # Dulcamaroid
  MRCA(tr4,"S._lyratum_Dulcamaroid", "S._aligerum_Dulcamaroid")
  # Morelloid
  MRCA(tr4,"S._radicans_Morelloid_Radicans", "S._zuloagae_Morelloid_Sect_Solanum")
  # Archaesolanum
  MRCA(tr4,"S._symonii_Archaesolanum", "S._laciniatum_Archaesolanum")
  # Normania
  MRCA(tr4,"S._trisectum_Normania", "S._herculeum_Normania")
  # African non-spiny
  MRCA(tr4,"S._imamense_African non???spiny", "S._guineense_African non???spiny")
  # Valdiviense
  MRCA(tr4,"S._valdiviense_Valdiviense", "S._alphonsei_Dulcamaroid")  
  
#Potato (This is with Regmandra)
MRCA(tr4,"S._piurae_Petota", "S._paposanum_Regmandra")

  #Regmandra
  #Pteroideae
  MRCA(tr4,"S._multifidum_Regmandra", "S._remyanum_Regmandra")
  #Herpystichum
  MRCA(tr4,"S._limoncochaense_Herpystichum", "S._loxophyllum_Herpystichumm")
  #Oxycoccoides
  #MRCA(tr4,"S._oxycoccoides_Oxycoccoides", "S._brevifolium_Anarrhichomenum")
  #Articulatum
  MRCA(tr4,"S._taeniotrichum", "S._sanctaemarthae_Articulatum")
  #Basarthrum
  MRCA(tr4,"S._suaveolens_Basarthrum", "S._cochoae_Basarthrum")
  #Anarrhicomenum
  MRCA(tr4,"S._complectens_Anarrhichomenum", "S._brevifolium_Anarrhichomenum")
  #Etuberosum
  MRCA(tr4,"S._palustre_Etuberosum", "S._etuberosum_Etuberosum")
  #Tomato
  MRCA(tr4,"S._ochranthum_Tomato", "S._chilense_Tomato")
  #Petota
  MRCA(tr4,"S._piurae_Petota", "S._oxycarpum_Petota")
  

#Clandestinum-Mapiriense
  #Anomalotemon
  #S._anomalostemon_Mapiriense
  #S._clandestinum_Mapiriense, S._mapiriense_Mapiriense
  
#Wendlandii-Allophyllum
  #S. graveolens
  #S._graveolens_Wendlandii???Allophyllum
  
  #S. allophyllum
  #Wendlandii-Allophyllum s.s.
  MRCA(tr4,"S._allophyllum_Wendlandii???Allophyllum", "S._bicorne_Wendlandii???Allophyllum")#860
  
  
#Nemorense
MRCA(tr4,"S._reptans_Nemorense", "S._nemorense_Nemorense")
  
#Cyphomandra
MRCA(tr4,"S._diploconos_Pachyphylla", "S._cacosmum_Pachyphylla")
  
  #Cyphomandropsis core
  MRCA(tr4,"S._amotapense_Cyphomandropsis", "S._stuckertii_Cyphomandropsis")

  #Pachylla
  MRCA(tr4,"S._pendulum_Pachyphylla", "S._cacosmum_Pachyphylla")
  MRCA(tr4,"S._cylindricum_Pachyphylla", "S._diploconos_Pachyphylla")
  
  #S._morellifolium_Cyphomandropsis
  #S._glaucophyllum_Cyphomandropsis

#Geminata
MRCA(tr4,"S._delitescens_Reductum", "S._acuminatum_Geminata")#860
  #Reductum
  MRCA(tr4,"S._delitescens_Reductum", "S._reductum_Reductum")#860
  #Geminata s.s.
  MRCA(tr4,"S._havanense_Geminata", "S._acuminatum_Geminata")#860
  
  
#Brevantherum
MRCA(tr4,"S._apiahyense_Brevantherum s.s.", "S._goodspeedii_Brevantherum s.s.")#860
  #Inornatum
  MRCA(tr4,"S._inornatum_Inornatum", "S._bradei_Inornatum")#860
  #Gonatotrichum
  MRCA(tr4,"S._turneroides_Gonatotrichum.", "S._lignescens_Gonatotrichum")#860
  
#Leptostemonum
MRCA(tr4,"S._polygamum_Leptostemonum_unplaced", "S._sessiliflorum_Lasiocarpa")#860

####
#Label Minor Clades

#CAMPECHIENSE
#S._campechiense_Campechiense

# CROTONOIDEs
MRCA(tr4,"S._crotonoides_Crotonoides", "S._woodburyi_Crotonoides")#860

# Old World
MRCA(tr4,"S._perplexum_Carolinense", "S._aridum_Carolinense")#860

# Carolinense
MRCA(tr4,"S._perplexum_Carolinense", "S._aridum_Carolinense")#860

# Hieronymi
#S._hieronymi_Hieronymi

# Asterophorum
MRCA(tr4,"S._asterophorum_Asterophorum", "S._piluliferum_Asterophorum")#860

# Bahamense
MRCA(tr4,"S._bahamense_Bahamense", "S._ensifolium_Bahamense")#860

# Multispinum
#S._multispinum_Multispinum

# Crinitum
MRCA(tr4,"S._mitlense_Crinitum", "S._urticans_Crinitum")#860

# Androceras
MRCA(tr4,"._cordicitum_Androceras", "S._rostratum_Androceras")#860

# Sisymbriifolium
MRCA(tr4,"S._hasslerianum_Sisymbriifolium", "S._sisymbriifolium_Sisymbriifolium")#860
  # S._euacanthum_Sisymbriifolium
  # S._reineckii_Sisymbriifolium

# Lasiocarpa
MRCA(tr4,"S._repandum_Lasiocarpa", "S._stramoniifolium_Lasiocarpa")#860

# Acanthophora
MRCA(tr4,"S._mammosum_Acanthophora", "S._affine_Acanthophora")#860

# Gardneri
MRCA(tr4,"S._stenandrum_Gardneri", "S._tetramerum_Gardneri")#860

# Thomasiifolium
MRCA(tr4,"S._paraibanum_Thomasiifolium", "S._buddleifolium_Thomasiifolium")#860

# Erythrotrichum
MRCA(tr4,"S._robustum_Erythrotrichum", "S._paludosum_Erythrotrichum")#860

# Torva
MRCA(tr4,"S._donianum_Torva", "S._pseudosaponaceum_Torva")#860

# Micracantha
MRCA(tr4,"S._jamaicense_Micracantha", "S._volubile_Micracantha")#860

# Elaeagniifolium
MRCA(tr4,"S._homalospermum_Elaeagnifolium", "S._houstonii_Elaeagnifolium")#860


###

#write posterior probability vlaues


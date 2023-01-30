#Evolution of the Neustonic Ecosystem
#Character Trait Reconstruction
#Anthony CJ, Bentlage B, Helm RR

#Warning: Names are often overwritten to avoid excessive objects in environment

#Install and load necessary packages
install.packages("readr")
library(readr)

install.packages("phylobase")
library(phylobase)
install.packages("treeio")
library(treeio)
install.packages("ape")
library(ape)
install.packages("phytools")
library(phytools)
install.packages("TreeTools")
library(TreeTools)
install.packages("phytools")
library(phytools)
install.packages("adephylo")
library(adephylo)

install.packages("ggplot2")
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("ggnewscale")
library(ggnewscale)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("gridExtra")
library(gridExtra)

#Load in RAxML tree
SiphoBoots <- read.raxml("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/RAxML_bipartitionsBranchLabels.Sipho_16S18S28S")
summary(SiphoBoots)

#Bootstrap value is visualized as color
#Node label is used to reroot tree
ggtree(SiphoBoots, aes(color=bootstrap)) + theme_tree2() +
  scale_color_continuous(high='#D55E00', low='#0072B2') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=node), hjust=-.3)

#Node 32 is our reroot
RootedTree <- root(SiphoBoots, node = 32, resolve.root = TRUE)

#Visualize rerooted tree
B <- ggtree(RootedTree, aes()) + theme_tree2() +
  geom_tiplab(size = 2)+
  geom_text(aes(label=bootstrap/100), hjust=-.3)+
  theme(legend.position = "none")
B

#Visualize ancestral states for categorical trait
#Use .csv with the exact trait and name

SiphoTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Siphonophore.csv", row.names=1)

#Remove any empty rows. They sometimes get loaded in at the bottom of the spreadsheet
#This code also renames our dataframe as 'trait_data'
SiphoTrait <- na.omit(SiphoTrait)

#change .csv to matrix fromat
head(SiphoTrait)
x<-as.matrix(SiphoTrait)[,1]
x

#change to phylo type tree
class(RootedTree)
phylotree <- as.phylo(RootedTree)
phylotree

cols<-setNames(c("black","#009e73ff","#d81b60ff","#1e88e5ff"),levels(x))

fitARD<-ace(x, phylotree, method = "ML", model="ER",type="discrete", marginal = TRUE)
fitARD

#ML values
fitARD$lik.anc
#Each character has a likelihood assigned to each node.
#Filtering through these values in conjunction with our original ggtree on row 38 will get us exact values

#Let's get an initial visualization of this reconstruction
plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD$lik.anc,piecol=cols,cex=1)

#While we could stop here, we can make this tree more informative and visually appealing
#Let's do that

#Reformat your tree as a "tibble"
xyztree <- as_tibble(RootedTree)
summary(xyztree)

#load in character trait data, this time making names their own column
SiphoTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Siphonophore.csv")
head(xyztree)
head(SiphoTrait)

#Use your orignial trat_data base to assign traits to taxa names
#This is where it is important they have the exact same name
z <- full_join(xyztree, SiphoTrait, by = 'label')
#Reformat this new 'y' object as .treedata
#This data form let's us maintain the character traits assigned to the tips
#Not all tree data types do this so be careful
z <- as.treedata(z)

#Make sure tip labels are correct
p<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Habitat), alpha = 0.6)+
  scale_color_manual(values = c("black","#009e73ff","#d81b60ff","#1e88e5ff"))
p

##Now we need to adjust your pies so they will map to object 'p'

#First let's make your ancestral states a dataframe
ancstats <- as.data.frame(fitARD$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)

#we only keep any pie charts that are informative
#In this case I kep anything with a probability of less than 0.99999 for each trait

ancstatssubset<- subset(ancstats, Pelagic <= c(0.99))

#We can get an idea of how many pie charts that removed with nrow
#This tells us how many rows are in our dataset

nrow(ancstats)
#Only 4 pies have likelihoods < 0.99 for pelagic trait
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(x))

#Then we assign colors to our character traits
#Make sure they are the same as your tip labels

#Sometimes they will get reversed or reorganized, so it may take some playing
pies <- lapply(pies, function(g) g + scale_fill_manual(values = c("black","#009e73ff","#d81b60ff","#1e88e5ff")))

#Now put your pies ontop of your already plotted tree (object 'p')
#width and height adjust the size of your pies, so adjust them to your liking
SiphoRecon <- p + geom_inset(pies, width = 0.05, height = 0.05)
SiphoRecon

#After you have a tree you like, I suggest exporting it as a .pdf
#This will maintain all of your vectors and keep the tree fully editable in post
#Typically I will take the tree from here, and import it into Inkscape
#There, I can remove uninformative pies, add informative ones, fix names, annotate, etc.

##
#FROM THIS POINT FORWARD, IT IS THE ABOVE ANNOTATED PIPELINE REPEATED FOR EACH GROUP
#NOTES WILL BE LIMITED

##
##Barnacles
##

#Load tree
BarnacleBoots <- read.raxml("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/RAxML_bipartitionsBranchLabels.Pedunc_18S28SCOI")

#Visualize
ggtree(BarnacleBoots, aes(color=bootstrap)) + theme_tree2() +
  scale_color_continuous(high='#000000ff', low='#ffffffff') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=node), hjust=-.3)

#Reroot
RootedTree <- root(BarnacleBoots, node = 41, resolve.root = TRUE)

#Check new root
A <- ggtree(RootedTree, aes()) + theme_tree2() +
  scale_color_continuous(high='#000000ff', low='#ffffffff') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=bootstrap/100), hjust=-.3)
A

#Load character traits
BarnacleTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Barnacles.csv", row.names=1)
head(BarnacleTrait)

x<-as.matrix(BarnacleTrait)[,1]
x

phylotree <- as.phylo(RootedTree)
phylotree

#Set colors
cols<-setNames(c("black","#009e73ff", "#d81b60ff", "#ffc107ff"),levels(x))

#Build model
fitARD<-ace(x, phylotree, method = "ML", model="ER",type="discrete", marginal = TRUE)
fitARD

fitARD$lik.anc

#Initial visualization of reconstruction
plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD$lik.anc,piecol=cols,cex=1)

#Beautification
xyztree <- as_tibble(RootedTree)
summary(xyztree)

#Load in character trait data with column label
BarnacleTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Barnacles.csv")
head(xyztree)
head(BarnacleTrait)

#Join tree and traits
z <- full_join(xyztree, BarnacleTrait, by = 'label')
z <-as.treedata(z)

#Make sure tip labels are correct
p<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Habitat), alpha = 0.6)+
  scale_color_manual(values = cols)
p

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)

#Keep informative pies
ancstatssubset<- subset(ancstats, Epibiotic <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Benthic <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Rafting <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(x))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols))
BarnacleRecon <- p + geom_inset(pies, width = 0.05, height = 0.05)
BarnacleRecon

##
##Capitata
##

CapiBoots <- read.raxml("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/RAxML_bipartitionsBranchLabels.Capitata_16S18S28S_")

ggtree(CapiBoots, aes(color=bootstrap)) + theme_tree2() +
  scale_color_continuous(high='#D55E00', low='#0072B2') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=node), hjust=-.3)

RootedTree <- root(CapiBoots, node = 62, resolve.root = TRUE)

C <- ggtree(RootedTree, aes()) + theme_tree2() +
  scale_color_continuous(high='#000000ff', low='#ffffffff') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=bootstrap/100), hjust=-.3)+
  theme(legend.position = "none")
C

CapiTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Capitata.csv", row.names=1)
head(CapiTrait)

#Make character matrix for both traits: Habitat and LifeCycle
x<-as.matrix(CapiTrait)[,1]
y<-as.matrix(CapiTrait)[,2]

phylotree <- as.phylo(RootedTree)

cols1<-setNames(c("black","#009e73ff","#d81b60ff","#1e88e5ff", "#b2b2b2ff"),levels(x))

fitARD<-ace(x, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD

fitARD$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD$lik.anc,piecol=cols,cex=0.5)

#Life Cycle
cols2<-setNames(c("black", "#d55e00ff", "#0e88e5ff"),levels(y))
fitARD2<-ace(y, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD2

fitARD2$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD$lik.anc,piecol=cols2,cex=0.5)

#Beautification
xyztree <- as_tibble(RootedTree)
summary(xyztree)

#Load in character trait data with column label
CapitataTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Capitata.csv")
head(xyztree)
head(CapitataTrait)

#Join tree and traits
z <- full_join(xyztree, CapitataTrait, by = 'label')
z <-as.treedata(z)

#Make sure tip labels are correct
p1<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Habitat), alpha = 0.6)+
  scale_color_manual(values = cols1)
p1

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Epibiotic <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Benthic <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(x))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols1))
CapitateRecon1 <- p1 + geom_inset(pies, width = 0.05, height = 0.05)
CapitateRecon1

#Make sure tip labels are correct
p2<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = LifeCycle), alpha = 0.6)+
  scale_color_manual(values = c("black", "#0e88e5ff", "#d55e00ff"))
p2

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD2$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Meroplanktonic <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Benthic <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(y))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = c("black", "#0e88e5ff", "#d55e00ff")))
CapitateRecon2 <- p2 + geom_inset(pies, width = 0.05, height = 0.05)
CapitateRecon2

##
##Cladobranchs
##

#Load in tree
CladoBoots <- read.raxml("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/RAxML_bipartitionsBranchLabels.CladoConstrained")
print(CladoBoots, n = 49)

#Plot tree with node labels
ggtree(CladoBoots, aes(color=bootstrap)) + theme_tree2() +
  scale_color_continuous(high='#D55E00', low='#0072B2') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=node), hjust=-.3)

#Reroot
RootedTree <- root(CladoBoots, node = 84, resolve.root = TRUE)

#Plot rerooted tree
D <- ggtree(RootedTree, aes()) + theme_tree2() +
  scale_color_continuous(high='#000000ff', low='#ffffffff') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=bootstrap/100), hjust=-.3)+
  theme(legend.position = "none")
D

#Load trait data
CladoTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Cladobranchia.csv", row.names=1)
head(CladoTrait)

#Make two seperate matrices
x<-as.matrix(CladoTrait)[,1]
x
y<-as.matrix(CladoTrait)[,2]
y

#Change object class
phylotree <- as.phylo(RootedTree)
phylotree

cols1<-setNames(c("black","#d81b60ff","#ffc207ff"),levels(x))

#Habitat
fitARD1<-ace(x, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD1

fitARD1$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD1$lik.anc,piecol=cols,cex=1)

#Prey
cols2<-setNames(c("#d671a4ff","#009e73ff","#d55e00ff","black", "#ffc207ff","#56b4e9ff"),levels(y))

fitARD2<-ace(y, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD2

fitARD2$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD2$lik.anc,piecol=cols2,cex=1)

#Beautification
xyztree <- as_tibble(RootedTree)
summary(xyztree)

#Load in character trait data with column label
CladobranchiaTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Cladobranchia.csv")
head(xyztree)
head(CapitataTrait)

#Join tree and traits
z <- full_join(xyztree, CladobranchiaTrait, by = 'label')
z <-as.treedata(z)

#Make sure tip labels are correct
p1<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Habitat), alpha = 0.6)+
  scale_color_manual(values = cols1)
p1

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD1$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Benthic <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(x))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols1))
CladobranchiaRecon1 <- p1 + geom_inset(pies, width = 0.05, height = 0.05)
CladobranchiaRecon1

#Make sure tip labels are correct
p2<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Prey), alpha = 0.6)+
  scale_color_manual(values = c(cols2))
p2

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD2$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Hydrozoa <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Hexacorallia <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(y))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols2))
CladobranchiaRecon2 <- p2 + geom_inset(pies, width = 0.05, height = 0.05)
CladobranchiaRecon2

##
##Epitoniids
##

EpiBoots <- read.raxml("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/RAxML_bipartitionsBranchLabels.Epi_16S28SH3")
print(EpiBoots, n = 15)
ggtree(EpiBoots, aes(color=bootstrap)) + theme_tree2() +
  scale_color_continuous(high='#D55E00', low='#0072B2') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=node), hjust=-.3)

RootedTree <- root(EpiBoots, node = 27, resolve.root = TRUE)

E <-ggtree(RootedTree, aes()) + theme_tree2() +
  scale_color_continuous(high='#000000ff', low='#ffffffff') + 
  geom_tiplab(size=2)+
  geom_text(aes(label=bootstrap/100), hjust=-.3)+
  theme(legend.position = "right")
E
EpiTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Epitoniidae.csv", row.names=1)
head(EpiTrait)

x<-as.matrix(EpiTrait)[,1]
x
y<-as.matrix(EpiTrait)[,2]
y

#Change object class
phylotree <- as.phylo(RootedTree)
phylotree

cols1<-setNames(c("black","#d81b60ff"),levels(x))

#Habitat
fitARD1<-ace(x, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD1

fitARD1$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD1$lik.anc,piecol=cols,cex=1)

#Prey
cols2<-setNames(c("#d55e00ff","black","#0072b2ff","#b2b2b2ff"),levels(y))

fitARD2<-ace(y, phylotree, method = "ML", type="discrete", marginal = TRUE)
fitARD2

fitARD2$lik.anc

plotTree(phylotree)
nodelabels(node=1:phylotree$Nnode+Ntip(phylotree),
           pie=fitARD2$lik.anc,piecol=cols2,cex=1)

#Beautification
xyztree <- as_tibble(RootedTree)
summary(xyztree)

#Load in character trait data with column label
EpitoniidaeTrait <- read.csv("/Users/colin/OneDrive/Documents/Neuston/NeustonCharacterTraitReconstruction-20230129T222802Z-001/NeustonCharacterTraitReconstruction/PieCharts/Epitoniidae.csv")
head(xyztree)
head(EpitoniidaeTrait)

#Join tree and traits
z <- full_join(xyztree,EpitoniidaeTrait, by = 'label')
z <-as.treedata(z)

#Make sure tip labels are correct
p1<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Habitat), alpha = 0.6)+
  scale_color_manual(values = cols1)
p1

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD1$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Benthic <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(x))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols1))
EpitoniidaeRecon1 <- p1 + geom_inset(pies, width = 0.05, height = 0.05)
EpitoniidaeRecon1

#Make sure tip labels are correct
p2<-ggtree(z) + theme_tree2()+
  geom_tiplab(size = 1)+
  theme(legend.position = "none")+
  geom_tippoint(aes(color = Prey), alpha = 0.6)+
  scale_color_manual(values = c(cols2))
p2

#Map pies to object 'p'
ancstats <- as.data.frame(fitARD2$lik.anc)
ancstats$node <- 1:phylotree$Nnode+Ntip(phylotree)
head(ancstats)

#Keep informative pies
ancstatssubset<- subset(ancstats, Hydrozoa <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Actiniidae <= c(0.99))
ancstatssubset<- subset(ancstatssubset, Scleractinia <= c(0.99))

nrow(ancstats)
nrow(ancstatssubset)

#Now we send our subset pie database to object "pies"
pies <- nodepie(ancstatssubset, cols = unique(y))
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols2))
EpitoniidaeRecon2 <- p2 + geom_inset(pies, width = 0.05, height = 0.05)
EpitoniidaeRecon2

###FINAL VISUALIZATIONS
grid.arrange(A,B,C,D,E)
grid.arrange(BarnacleRecon, SiphoRecon)
grid.arrange(CapitateRecon1, CapitateRecon2)
grid.arrange(CladobranchiaRecon1, CladobranchiaRecon2)
grid.arrange(EpitoniidaeRecon1, EpitoniidaeRecon2)

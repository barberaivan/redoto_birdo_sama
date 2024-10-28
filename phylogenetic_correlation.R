#############################################
##phylogenetic distances
#birds
library(ape)

treebird<- read.nexus("birds_tree.nex") #load the phylogenetic tree

plot.phylo(treebird[[55]], type= "fan") #let's see how horrible it is

abird<- multi2di(treebird)#solve the mutichotomies with a random procedure

##take the phylogenetic correlation matrices for all the hypothetic trees
treelist = vector('list', 3)
for (i in 1:length(abird)){
  treelist[[i]] = as.matrix(vcv.phylo(abird[[i]], cor=TRUE))
}
## average the phylogenetic correlation into one matrix
bircor = apply(simplify2array(treelist), 1:2, mean)


##plants

plantas<-read.table("plantas.txt", sep=";", header = TRUE)
plant<- as.data.frame(unique(plantas$plant.sp))
colnames(plant)<- "plant.sp"
genus<- matrix(nrow= nrow(plant))
for( i in 1:nrow(plant)){
  tmp = strsplit(as.character(plant$plant.sp[i]), " ") ##take the genus
  genus[i,] = sort(substring(tmp[[1]][1], 1, ))

}


devtools::install_github("ropensci/taxize")##need this to get the families
library(taxize)
plant$genus<- genus


plant$family <- tax_name(plant$genus, get = 'family', db = 'ncbi')$family  ##get the families


devtools::install_github("jinyizju/V.PhyloMaker2") ###get the megatree of plants and the tools to prune it

library(V.PhyloMaker2)

plantree<- phylo.maker(plant, scenarios = "S1")##use the plant data frame to get the trees
plot.phylo(plantree[[1]], type="fan")#see how horrible it is
plantree<- multi2di(plantree[[1]])#solve multichotomies (luckilly they are not so much)

placor<- as.matrix(vcv.phylo(plantree, cor=TRUE)) #get the phylogenetic correlation matrix

save(bircor, placor, file="filogenias.R") #save both birds and plant matrices
load("filogenias.R") #load them to happyly analyse the data
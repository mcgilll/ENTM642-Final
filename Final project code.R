
######  CAWS DATA 5 sites per waterway  ######
install.packages("vegan")
install.packages("ape")
install.packages("dendextend")
library(vegan)									
library(ape)									### PCoA package
library(dendextend)								### Dendrogram package

###### Read in and clean up data ####

CCDS<-read.csv("CAWS fish abundance data.csv")
str(CCDS)
CCDS[is.na(CCDS)] <- 0								### Replaced NA values with zeros
str(CCDS)										### Checked to make sure
dim(CCDS)
my.means <- colMeans(CCDS[,3:79])
my.means										### Wide range of species abundance
str(my.means)
common.species <- subset(CCDS[3:79], select=which(my.means>1))
colMeans(common.species)							### Subset df to only include common species

scaled.ccd<-sqrt(common.species)						### I scaled down the data because some species are much more abundant than others. 
									
scaled.ccd										### Square root of original data
	
######  PCOA  #######
waterways <-as.factor(CCDS[,1])							### Made waterways a factor
caws.full <- cbind(scaled.ccd, waterways)						### Combined waterways and scaled fish data
dim(caws.full)


fish.d.bc <- vegdist(caws.full[,1:37], "bray")					### Created dissimilarity matrix using bray-crtis dissimilarity
fish.bc.pcoa <- pcoa(fish.d.bc)							### Created PCoA from dissimilarity matrix
biplot(fish.bc.pcoa)									### Plot PCoA

fish.d.j <- vegdist(scaled.ccd, "jaccard")					### Created dissimilarity matrix using jaccard dissimilarity based on presence/absence
fish.pcoa.j <- pcoa(fish.d.j)
biplot(fish.pcoa.j)

par(mfrow=c(1,2), cex=0.7)
biplot(fish.bc.pcoa, main = "Bray-Curtis PCoA ordination")
biplot(fish.pcoa.j, main= "Jaccard PCoA ordination")
########## Dendrogram   ########
									

clust.j.w <- hclust(fish.d.j, method="ward.D2")							##  Use 2 different cluster methods 
clust.bc.w <- hclust(fish.d.bc, method="ward.D2")

clust.bc.w <- as.dendrogram(clust.bc.w)
clust.j.w <-as.dendrogram(clust.j.w)
par(mfrow=c(2,1),cex=0.6, col="darkblue", mar=c(6,2,2,10))						### Plotting dendrogram of both jaccard and B-C to compare
plot(clust.bc.w, horiz=TRUE, main="Bray Curtis", xlab="Distance in Species Space")		### Very similar
plot(clust.j.w, horiz=TRUE, main="Jaccard", xlab="Distance in Species Space")
	
par(mfrow=c(1,1), cex=1)

#####  Testing for differences between waterways ####

waterway <- as.factor(CCDS[,1])								## Separated the waterways. Made them a factor
fish.full <- cbind(scaled.ccd, waterway) 							## Combined species data and waterways into one dataframe
attach(fish.full)
str(fish.full)
species <- fish.full[,1:37]									## Separated the species					

adon.disp <- betadisper(vegdist(fish.full[,1:37], method="bray"), waterway)	## Check assumption of similar dispersions
anova(adon.disp)											## Group dispersions are homogenous

waterways.adon <- adonis(species ~ waterway, method="bray", data=fish.full,	## Do adonis to check for differences
	control=permControl(strata=waterway), permuations=9999)		 	## Test against 9999 permutations as null
waterways.adon											## Check results. They are significant
plot(adon.disp)											## Plot the dispersions. NSC is significantly different from SBCR and NBCR
	
########################



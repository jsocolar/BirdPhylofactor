##Phylofactorization of Costa Rican birds data

### If you haven't already installed the updated phylofactor package, run the 5 lines below
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# biocLite("Biostrings")
# install.packages('devtools')
# devtools::install_github('reptalex/phylofactor')

library(phylofactor)

### import & pre-process data
X <- read.csv('data/CRBirds.csv',header=T,
              colClasses = c('factor','factor','numeric','factor','character','numeric')) %>%
  as.data.table
X <- X[!duplicated(X),]   #removing duplicate rows - not sure what to do with those.
names(X)[names(X)=='Spp'] <- 'Species'
X[,Species:=gsub(' ','_',Species)]
bird.species <- unique(X$Species)
X[,LUI:=factor(LUI,levels=c('Forest','Diversified','Intensive'))]

taxonomy <- read.csv('data/2012-03-04206D-master_taxonomy.csv',header=T,stringsAsFactors = F)
tiplabels <- taxonomy$TipLabel
tax <- taxonomy[,c('OscSubOsc','Order','BLFamilyLatin','Clade','TipLabel')]
tax <- apply(tax,1,paste,collapse='; ')   ##semicolon delimited taxonomies are necessary for taxonomic summaries
taxonomy <- data.frame("Species"=tiplabels,'taxonomy'=tax,stringsAsFactors = F)

### Start with a single tree
tree <- read.tree('data/JETZ TREES/Jetz_with_Furnariidae_March2013__tree1.txt')


bird.species[!bird.species %in% tree$tip.label] #18 species not found in tree tip labels
# [1] "Veniliornis_fumigatus"     "Pipra_coronata"           
# [3] "Xenops_minutus"            "Contopus_sp."             
# [5] "Pionopsitta_haematotis"    "Piculus_rubiginosus"      
# [7] "Automolus_rubiginosus"     "Buarremon_brunneinucha"   
# [9] "Myrmotherula_fulviventris" "Ciccaba_virgata"          
# [11] "Amazilia_sp."              "Empidonax_sp."            
# [13] "Laterallus_sp."            "Ardea_alba"               
# [15] "Ceryle_torquatus"          "Gallinago_delicata"       
# [17] "Ceryle_alcyon"             "Buarremon_torquatus" 

X <- X[Species %in% tree$tip.label,]  ## filter dataset to species in tree
bird.species <- unique(X$Species)
tree <- drop.tip(tree,setdiff(tree$tip.label,bird.species)) ## trim tree to only species in data


X[,Sample:=paste(Transect,Year,Season,LUI,sep='_')]
X[,N:=Number.of.surveys.per.season.detection..out.of.3.]

### There are two main for generalized phylofactorization of these data
## (1) algorithm='mStable' - sum counts within groups and use factor contrasts to ID edges best separating species
## (2) algorithm='phylo' - skip within-group summation. Much more costly

### We'll use the mStable algorithm.
## While the function gpf() accepts glm-style input of data frames containing all the variables
## used in regression, the mStable algorithm will convert that data frame into a matrix & operate on it
## accordingly. To make our workspace easily merged with the mStable phylofactor object,
## we'll make those matrices. For binoimal regression, need a list containing two matrices: Successes & Failures.
## In this script, the absence of a species is considered 0 successes, 3 failures.

Successes <- phyloframe.to.matrix(X)
## We'll put the successes & failures in the same order as tree$tip.labels
Successes <- Successes[tree$tip.label,]
Failures <- 3-Successes

### The meta-data must have rows match each column in Successes.
MetaData <- X[,c('Sample','Season','LUI','Year','Transect')]
MetaData <- MetaData[!duplicated(MetaData)]
MetaData <- MetaData[match(colnames(Successes),Sample),]


Data <- list('Successes'=Successes,'Failures'=Failures)
pf <- gpf(Data,tree,MetaData=MetaData,frmla.phylo=cbind(Successes,Failures)~Season+phylo*LUI,
          nfactors=200,ncores=7,algorithm = 'mStable',family=binomial,PartitioningVariables = 'LUI')
## If we don't set PartitioningVariables='LUI', then gpf will also partition edges based on the deviance
## from the term phylo all by itself (i.e. edges separating species with a high difference in P{Success} 
## irrespective of how that P{Success} changes with LUI)

MetaData[,LUI_quant:=as.numeric(LUI)]

## We can also use a quantitative model of increasing land-use intensity. The benefit being
## that, while the phylo*LUI factor above splits differences-of-means in each group, the quant model
## focuses more on the response slopes to LUI, not the means (i.e. relative difference-of-means)

pf_quant <- gpf(Data,tree,MetaData=MetaData,frmla.phylo=cbind(Successes,Failures)~Season+phylo*LUI_quant,
              nfactors=200,ncores=7,algorithm = 'mStable',family=binomial,PartitioningVariables = 'LUI_quant')

save(list=ls(),file='workspaces/CR_LUI_mStable_phylofactorization_workspace')

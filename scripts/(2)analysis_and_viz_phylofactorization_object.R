library(phylofactor)
library(ggplot)
library(ggpubr)

load('workspaces/CR_LUI_mStable_phylofactorization_workspace')




# Summary -----------------------------------------------------------------

pf
### the following function will be useful for visualizing results from each factor

pf.boxplot <- function(n,pf.=pf){
  pf$models[[n]]$data %>%
      ggplot(aes(LUI,Successes,color=phylo))+
      geom_boxplot()+
      geom_jitter(size=3,width = .2)+
      facet_grid(.~phylo)+
      ggtitle(paste('Factor',n))
}

### an easy way to peruse phylofactor objects
n=7
s <- summary(pf,taxonomy,n)                                     ### shortest-unique-taxonomic-prefix grouping & signal calculation for relative importance
s.species <- summary(pf,taxonomy,n,taxon.trimming='species')  ### provides species-specific signal calculation

s$taxa.split$group1 ##group1 is the descendant and often monophyletic group, R, in the data below:
s.species
#### The species-level summary re-computes the objective function for each species, contrasting it
#### with the entirety of the opposing group as a measure of "signal" for each species. Thus, if we call
#### Empidonax_flaviventris R, all of Group2 S, and all the rest of species NA, we get 
#### a deviance of 742.0165 from the phylo:LUI terms. Only the top-three species for each level are displayed
#### when printing the summary object, but the full information can be obtained in the signal.table element.
s.species$signal.table$Group1[1:10,]
s.species$signal.table$Group2[1:10,]

#### One can usually expect clades at the uppper & lower bounds of these signal tables to be
#### partitioned in later factors.
# For example, factors 8-10 further partition the paraphyletic Group2
# and factors 17 and 18 further partition the suboscine group in 7


pf$tree$tip.label[sapply(pf$groups[8:10],'[[',1)] 
#factors 8=10: Turdus grayi, Hylophilus_decurtatus, Henicorhina_leucosticta, top members of the Group2 signal.table

summary(pf,taxonomy,17)  # Synallaxis_brachyura
pf.boxplot(17)
summary(pf,taxonomy,18)
pf.boxplot(18)           #A group of TitTyranRest, led by Ornithion_brunneicapillus, don't have as dramatic of a decrease from forest-diversified




# Number of sign factors --------------------------------------------------

### Quick check: What is an approximate number of significant factors?
### To be conservative, we'll use a strict, bonferroni cutoff of 
### multiple comparisons, where the number of comparisons is the number of
### edges in the tree
n_edge <- Nedge(tree)

model_deviance <- numeric(pf$nfactors)
for (i in 1:pf$nfactors){
  model_deviance[i] <- pvDeviance(pf$models[[i]],pf$groups[[i]],pf$tree,
                                  PartitioningVariables = 'LUI',model.fcn = glm,pf$models[[i]]$data)
}
pvals <- 1-pchisq(model_deviance,df=2)
cutoff <- 0.05/n_edge   
min(which(pvals>cutoff)) 
#all of these factors are well below the cutoff.












# More Visualization -----------------------------------------------------------

tree_plot <-  pf.tree(pf,factors = 1:7)
tree_plot$ggplot

### Demo: how to plot data for a factor that matches colors in the tree_plot
n=7
box_plot <- pf$models[[n]]$data[pf$models[[n]]$data$phylo=='R'] %>%
  ggplot(aes(LUI,Successes))+
  geom_boxplot()+
  geom_jitter(color=tree_plot$legend$colors[n],size=3,width = .2)+
  ggtitle(paste('Factor',n))

ggarrange(tree_plot$ggplot,box_plot,ncol=2)

### for ease, we'll make a function that box-plots a particular factor
### given a pf.tree object
pp <- pf.tree(pf,factors = 1:10)
bplot <- function(i){
  if (i<=nrow(pp$legend)){  ## the paraphyletic remainder is stored in the last row
    grp <- pf$groups[[i]][[1]]
    col <- pp$legend$colors[i]
  } else {
    grp <- pf$groups[[i]][[2]]
    col <- NULL
  }
  y <- colSums(Successes[grp,,drop=F])
  plot(X$LUI,y,xlab=NULL,ylab=NULL,lwd=2,cex.axis=2,col=col,cex.main=1.5,main=paste('Factor',i))
}


tiff('Figures/All_Birds_Counts.tiff',height=400,width=800)
  par(mfrow=c(1,1))
  y <- colSums(Successes)    #ag
  plot(X$LUI,y,xlab=NULL,ylab=NULL,lwd=2,main="All Birds",col='grey',cex.axis=2,cex.main=2)
dev.off()

tiff('Figures/Bird_Factor_Counts.tiff',height=1400,width=800)
  par(mfrow=c(3,2))
  for (i in 1:6){
    bplot(i)
  }
dev.off()


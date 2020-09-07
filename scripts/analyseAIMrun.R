library(ggplot2)
library(ape)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source("./plotSpeciesTree.R") # sources the function to plot the species tree with arrows
trees <- "./../precooked_runs/summary.tree" # defines the species tree location
species_trees <- read.nexus(trees) # get all the tree topologies in the posterior

prior = 0.695/50 # define the prior probability of an indicator being active (here it's lambda/(nr indicators) )
posterior.threshold = 0.5 # migration has to have a higher posterior support to be plotted
BF.threshold = 0 # migration has to have a higher Bayes Factor to be plotted
forwards.arrows = T # if true, the arrows denote forward in time migration, if false, backwards in time migration

# compute total number of trees, only works if Minimal tree support for output was 0 in AIM ann
tot.trees = 0
for (i in seq(1, length(species_trees))){
  post.occurance = strsplit(names(species_trees)[[i]], split="_")[[1]]
  tot.trees = tot.trees +as.numeric(post.occurance[[length(post.occurance)]])
}

for (i in seq(1, length(species_trees))){
  plot_tree <- ladderize(species_trees[[i]], right=F) # ladderizes the tree (defines the node orderign on the x-axis for the plot)
  p <- plotSpeciesTree(plot_tree, posterior.threshold, BF.threshold, prior, forwards.arrows) # returns a ggplot object that can be modified
  post.occurance = strsplit(names(species_trees)[[i]], split="_")[[1]] # get how often this ranked topology was observed in the posterior
  p <- p + 
    ylab("substitutions per site and year") + # labels the y-axis
    scale_size_continuous(range=c(0.1,1)) + # defines the range of arrow sizes
    ggtitle(paste("posterior support for ranked topology was:\n", as.numeric(post.occurance[[length(post.occurance)]])/tot.trees, sep="")) + # specify the title
    theme_minimal() # changes the plotting theme
  plot(p) # plots the species tree
  ggsave(paste(names(species_trees)[[i]], ".pdf", sep="")) 
}

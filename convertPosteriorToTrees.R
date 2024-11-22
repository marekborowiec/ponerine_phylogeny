library(stringr)
library(dplyr)
library(tidyr)
library(ape)
library(phangorn)
library(MCMCtreeR)

# Read in posterior file with node ages
MCMC_node_ages <- read.csv("./renumbered-concat-mcmc.txt", sep="\t")
# Read in the topology
MCMC_phy <- readMCMCtree("./FigTree.tre")
phy <- MCMC_phy$apePhy

# Source function
source("./functions/build_posterior_phylogenies_from_MCMC_outputs.R")

# Build new posterior phylogenies
PP_phylo <- build_posterior_phylogenies_from_MCMC_outputs(MCMC_node_ages = MCMC_node_ages, 
                                                      phy = MCMC_phy$apePhy,
                                                      time_conversion = 100)

# Set the number of trees to sample
num_trees_to_sample <- 1000

# Randomly select 1,000 tree indices
selected_indices <- sample(length(PP_phylo), num_trees_to_sample)

# Extract the selected trees
selected_trees <- PP_phylo[selected_indices]

# Define the output file name
output_file <- "1000_random_trees.txt"

# Write the selected trees to the file
write.tree(selected_trees, file = output_file)

# Function to print treePL config
write_treepl_config <- function(index, tree) {
    # Open connection to output file
    sink(paste("treepl.", "config", index, ".txt", sep=""))
    # Print input tree file name
    cat("treefile = rooted-no-support-full-tree.tre\n")
    #cat("smooth = 100\n")
    cat("numsites = 1055001\n")
    # Initialize a vector for node ages
    n_nodes <- Nnode(tree) + Ntip(tree)
    node_ages <- numeric(n_nodes)
    # Compute the age of each node
    for (i in 1:n_nodes) {
        age <- abs(round(node.depth.edgelength(tree)[i],5) - round(max(node.depth.edgelength(tree)),5))
        node_ages[i] <- age
        if (i > Ntip(tree)) {
            no_descendants <- length(tree$tip.label[Descendants(tree, i, type="tips")[[1]]])
            descendant1 <- tree$tip.label[Descendants(tree, i, type="tips")[[1]]][1]
            descendant2 <- tree$tip.label[Descendants(tree, i, type="tips")[[1]]][no_descendants]
            cat("mrca = node", i, " ", descendant1, " ", descendant2, "\n", sep="")
            cat("min = node", i, " ", age, "\n", sep="")
            cat("max = node", i, " ", age, "\n", sep="")
        }
    }
    # Write output file
    cat("outfile = tree", index, ".dated.tre\n", sep="")
    cat("thorough\n")
    cat("randomcv\n")
    cat("cvstart = 1\n")
    cat("cvstop = 0.0000001")
    # Close connection
    sink()
}

trees <- read.tree("1000_random_trees.txt")

for (tree_i in seq_along(trees)) {
    write_treepl_config(tree_i, trees[[tree_i]])
}
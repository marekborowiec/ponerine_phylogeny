library(ape)

# Load your phylogenetic tree from a newick-format file
tree <- read.tree("ponerinae-792t-spruce-75p-iqtree-swscmerge-mfp_v2_v2.tre")

# Strip branch lengths
tree$edge.length <- NULL

# Remove support values from all nodes
tree$node.label <- NULL

# Root
consensus_tree <- root(tree, c("Amblyopone_australis_D0872_CASENT0106229",
	"Paraponera_clavata_EX1573_CASENT0633292",
	"Proceratium_google_MAMI0434_CASENT0035028"))

mcmctree_annotator <- function(triplet, tree) {
	tip1 <- triplet[[1]]
	tip2 <- triplet[[2]]
	age <- triplet[[3]]
	# Find the MRCA node
	mrca_node <- getMRCA(tree, c(tip1, tip2))
	# Annotate the MRCA node with your text
	tree$node.label[mrca_node - length(tree$tip.label)] <- age
	return(tree)
}

calibrations <- list (
c("Odontoponera_denticulata_EX2968_CASENT0267032","Mesoponera_afr05_D2394_CASENT0815823",">.15"),
c("Neoponera_crenata_EX2451_CASENT0649900","Neoponera_foetida_EX2274_CASENT0633624",">.15"),
c("Platythyrea_arnoldi_D2118_CASENT0842276","Diacamma_geometricum_EX2682_ZRC_ENT00047849",">.94"),
c("Ponera_jtl001_EX2356_INB0003621410","Emeryopone_my01_D0990_CASENT0226984",">.34"),
c("Hypoponera_exigua_EX2769_CASENT0226547","Centromyrmex_hamulatus_EX2669_ZRC_ENT00047844",">.34"),
c("Odontomachus_hastatus_EX2290_CASENT0649099","Odontomachus_haematodus_EX2301_CASENT0649088",">.15"),
c("Anochetus_emarginatus_EX2286a_CASENT0646105","Anochetus_mayri_EX2271_CASENT0633416",">.15")
)

for (triplet in calibrations) {
	consensus_tree <- mcmctree_annotator(triplet, consensus_tree)
}

consensus_tree <- root(consensus_tree, c("Amblyopone_australis_D0872_CASENT0106229",
	"Paraponera_clavata_EX1573_CASENT0633292",
	"Proceratium_google_MAMI0434_CASENT0035028"), resolve.root=TRUE)

write.tree(consensus_tree, "MCMCTree_annotated_AoW.tre")

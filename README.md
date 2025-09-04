# Scripts for dating analyses

This repository contains custom scripts used for divergence dating analyses from the following publication:

`Doré M., Borowiec M.L., Branstetter M.G., Camacho G.P., Fisher B.L., Longino J.T., Ward P.S., & Blaimer B.B. (2025) Evolutionary history of ponerine ants highlights how the timing of dispersal events shapes modern biodiversity. Nature Communications. DOI: 10.1038/s41467-025-63709-3.`

### Scripts
`check_fossils.R` Script to verify that fossil calibrations are valid

`convertPosteriorToTrees.R` Wrapper script for build_posterior_phylogenies_from_MCMC_outputs.R to extract trees and print TreePL configuration files

`build_posterior_phylogenies_from_MCMC_outputs.R` Script to convert MCMCTree posterior files into tree files

`MCMCTree_creator_AoW.R` Script to automate creation of a tree file with node calibrations for MCMCTree input

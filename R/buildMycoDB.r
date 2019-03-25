#strains <- gsub("strain=", "", asm.sub3$infraspecific_name)
#tmp <- sapply(1:nrow(asm.sub3), function(i) {gsub(strains[i], "", asm.sub3$organism_name[i])})

#asm.sub3$species_name <- sapply(asm.sub3$organism_name,
#       function(organism_name) {
#         tmp <- sapply(speciesNames, grepl, organism_name)
#         if (any(tmp)) {
#           speciesNames[max(which(tmp))]
#         } else {
#           ""
#         }})

#o <- order(subset(asm.sub3, grepl("Mycobacterium abscessus", species_name))[,1])
#subset(asm.sub3, grepl("Mycobacterium abscessus", species_name))[o[1:5],1]


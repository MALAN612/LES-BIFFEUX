############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")

install.packages('ape')

install.packages("phangorn", dependencies = TRUE)


library(DECIPHER)
library(ape)
library(phangorn)
library(Biostrings)

########################################################
# ----- a) Importation et alignement (DECIPHER) ------ #
########################################################

# Importation des données
data <- readDNAStringSet("Cetacea COI.fasta")

# Alignement
aligned_data <- AlignSeqs(data)

# Retire les gaps de l'alignement pour retrouver les séquences brutes
sequences_brutes <- RemoveGaps(aligned_data)

# Alignement avec les pénalités à zéro
aligned_no_gap <- AlignSeqs(sequences_brutes, gapOpening = 0, gapExtension = 0)


# Visualiser l'alignement dans un browser (optionnel)
BrowseSeqs(aligned_no_gap, highlight = 0)

# Sauvegarde l'alignement final
writeXStringSet(aligned_no_gap, "aligned_sequences.fasta")

##################################################
# ----- b) Traduction et cadre de lecture ------ #
##################################################

# Traduction avec code mitochondrial vertébré
aligned_dna <- AlignTranslation(data,
                                readingFrame = 1,
                                geneticCode = getGeneticCode("2"),
                                type = "DNAStringSet")

# Pour visualiser l'alignement en AA
aligned_aa_1 <- AlignTranslation(data,
                                  readingFrame = 1,
                                  geneticCode = getGeneticCode("2"),
                                  type = "AAStringSet")
BrowseSeqs(aligned_aa_1)

# Test le cadre de lecture 2
aligned_aa_2 <- AlignTranslation(data,
                                  readingFrame = 2,
                                  geneticCode = getGeneticCode("2"),
                                  type = "AAStringSet")
BrowseSeqs(aligned_aa_2)

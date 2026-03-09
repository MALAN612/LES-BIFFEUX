############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

install.packages('ape')

install.packages("phangorn", dependencies = TRUE)


library(DECIPHER)
library(Biostrings)
library(ape)
library(phangorn)


########################################################
# ----- a) Importation et alignement (DECIPHER) ------ #
########################################################

# Importation des séquences
seqs <- readDNAStringSet("Cetacea COI.fasta")

# Alignement avec paramètres par défaut
alignment_default <- AlignSeqs(seqs)

# Alignement avec pénalités d’indels = 0
alignment_no_gap_penalty <- AlignSeqs(seqs, gapOpening = 0, gapExtension = 0)

# Export pour vérification
writeXStringSet(alignment_default, "alignment_default.fasta")
writeXStringSet(alignment_no_gap_penalty, "alignment_no_gap_penalty.fasta")

##################################################
# ----- b) Traduction et cadre de lecture ------ #
##################################################

# Alignement d'acides aminées
aa_alignment_1 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = 1)

writeXStringSet(aa_alignment, "aa_alignment_1.fasta")

aa_alignment_2 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = 2)

writeXStringSet(aa_alignment_2, "aa_alignment_2.fasta")

# Vérification
# Nombre de codons stop par séquence
stops_frame1 <- letterFrequency(aa_alignment_1, "*")

# Résumé statistique (min, médiane, max)
summary(stops_frame1)

# Distribution : combien de séquences ont 0, 1, 2 stops, etc.
table(stops_frame1)

# Identifier les séquences contenant au moins un codon stop
names(aa_alignment_1)[stops_frame1 > 0]

# Nombre de codons stop par séquence
stops_frame2 <- letterFrequency(aa_alignment_2, "*")

# Résumé statistique
summary(stops_frame2)

# Distribution du nombre de stops
table(stops_frame2)

# Identifier les séquences contenant au moins un codon stop
names(aa_alignment_2)[stops_frame2 > 0]

#############################################
# ------ c) Phylogénie nucléotidique ------ #
#############################################

# Alignement en format "DNAbin" pour ape
dna_bin <- as.DNAbin(alignment_default)

# Matrices de distance avec les bons codes
dist_JC   <- dist.dna(dna_bin, model = "JC69")
dist_K2P  <- dist.dna(dna_bin, model = "K80")
dist_TN93 <- dist.dna(dna_bin, model = "TN93")
dist_GG   <- dist.dna(dna_bin, model = "GG95") # Correction ici

# Force la conversion en matrice.
dna_mat <- as.matrix(dna_bin)

# Vérifie que c'est bien une matrice
is.matrix(dna_mat) # Doit retourner TRUE

# Analyse Bootstrap
run_bootstrap <- function(dna_input, mod) {
    # Arbre de référence
    d <- dist.dna(dna_input, model = mod)
    tree <- nj(d)

    # Bootstrap (1000 itérations)
    boot_scores <- boot.phylo(tree, dna_input, function(x) nj(dist.dna(x, model = mod)), B = 1000)

    return(boot_scores)
}

# Lance les calculs sur dna_mat
boot_JC   <- run_bootstrap(dna_mat, "JC69")
boot_K2P  <- run_bootstrap(dna_mat, "K80")
boot_TN93 <- run_bootstrap(dna_mat, "TN93")
boot_GG   <- run_bootstrap(dna_mat, "GG95")

# Calcul des moyennes en ignorant les NA
scores_moyens <- c(
    JC   = mean(boot_JC, na.rm = TRUE),
    K2P  = mean(boot_K2P, na.rm = TRUE),
    TN93 = mean(boot_TN93, na.rm = TRUE),
    GG   = mean(boot_GG, na.rm = TRUE)
)

print(scores_moyens)

# Identifie le modèle le plus robuste
meilleur_ape <- names(which.max(scores_moyens))
cat("Le modèle le plus robuste selon le bootstrap est :", meilleur_ape, "\n")


# Conversion au format phyDat
dna_pd <- as.phyDat(dna_mat)

# Test de modèles
mt <- modelTest(dna_pd)

# Trouve le meilleur modèle selon le critère BIC
best_model <- mt[which.min(mt$BIC), ]
print(best_model)

# Transforme le meilleur résultat du test en objet PML
fit_best <- as.pml(mt, model = "BIC")

alpha_val <- fit_best$shape
inv_val   <- fit_best$inv

cat("Paramètre Gamma (alpha) :", alpha_val, "\n")
cat("Sites invariants (I) :", inv_val, "\n")

# Optimisation de l'arbre
fit_optim <- optim.pml(fit_best, optGamma = TRUE, optInv = TRUE, optEdge = TRUE)

# Affichage de l'arbre final
par(mar = c(1, 1, 3, 1))
plot(fit_optim$tree, main = paste("Arbre Final (ML) - Modèle :", best_model$Model), cex = 0.8)
add.scale.bar()


###################################################################
# ----- d) Analyse des acides aminées et cadres de lecture ------ #
###################################################################

# Calcul pour l'alignement standard
dist_LG  <- dist.ml(aa_pd_default, model="LG")
dist_JTT <- dist.ml(aa_pd_default, model="JTT")
dist_B62 <- dist.ml(aa_pd_default, model="Blosum62")

# Fonction pour automatiser le bootstrap AA
run_aa_bootstrap <- function(pd_data, model_name) {
    # Arbre NJ de référence
    tree_ref <- nj(dist.ml(pd_data, model = model_name))

    # Bootstrap
    # On convertit en matrice pour boot.phylo
    boot <- boot.phylo(tree_ref, as.matrix(as.character(pd_data)),
                       function(x) nj(dist.ml(as.phyDat(x, type="AA"), model = model_name)),
                       B = 1000, quiet = TRUE)
    return(boot)
}

# Calcul de la robustesse pour le cadre par défaut (Standard)
boot_LG_std  <- run_aa_bootstrap(aa_pd_default, "LG")
boot_JTT_std <- run_aa_bootstrap(aa_pd_default, "JTT")
boot_B62_std <- run_aa_bootstrap(aa_pd_default, "Blosum62")

# Calcul pour le cadre de lecture décalé (LG est le témoin)
boot_LG_shift <- run_aa_bootstrap(aa_pd_shift, "LG")

# Calcul des moyennes en ignorant les NA (valeurs manquantes)
res_aa <- data.frame(
    Modèle = c("LG", "JTT", "Blosum62", "LG_Shifté"),
    Moyenne_Bootstrap = c(
        mean(boot_LG_std, na.rm = TRUE),
        mean(boot_JTT_std, na.rm = TRUE),
        mean(boot_B62_std, na.rm = TRUE),
        mean(boot_LG_shift, na.rm = TRUE)
    )
)

print(res_aa)

# Visualisation des deux topologies AA pour comparaison
par(mfrow = c(1, 2))
plot(nj(dist_LG), main = "Topologie Standard (LG)")
plot(nj(dist.ml(aa_pd_shift, model="LG")), main = "Topologie Shiftée (Bruit)")

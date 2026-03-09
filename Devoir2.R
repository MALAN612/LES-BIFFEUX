############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")

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

writeXStringSet(aa_alignment_1, "aa_alignment_1.fasta")

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

# Conversion de l'alignement en format DNAbin (ape)
dna <- as.DNAbin(alignment_default)

# Matrices de distances avec les différents modèles
dist_JC  <- dist.dna(dna, model = "JC69")    # Jukes-Cantor
dist_K80 <- dist.dna(dna, model = "K80")     # Kimura 2 paramètres
dist_TN  <- dist.dna(dna, model = "TN93")    # Tamura-Nei
dist_LOG <- dist.dna(dna, model = "logdet")  # Galtier-Gouy (LogDet)

# Construction des arbres
tree_JC  <- nj(dist_JC)
tree_K80 <- nj(dist_K80)
tree_TN  <- nj(dist_TN)
tree_LOG <- nj(dist_LOG)

# boot.phylo nécessite une matrice
dna_matrix <- as.matrix(dna)

# Bootstrap (1000 itérations)
boot_JC <- boot.phylo(tree_JC, dna_matrix,
                      function(x) nj(dist.dna(x, model="JC69")),
                      B = 1000)

boot_K80 <- boot.phylo(tree_K80, dna_matrix,
                       function(x) nj(dist.dna(x, model="K80")),
                       B = 1000)

boot_TN <- boot.phylo(tree_TN, dna_matrix,
                      function(x) nj(dist.dna(x, model="TN93")),
                      B = 1000)

boot_LOG <- boot.phylo(tree_LOG, dna_matrix,
                       function(x) nj(dist.dna(x, model="logdet")),
                       B = 1000)

# Calculer la moyenne en ignorant les valeurs manquantes
mean(boot_JC, na.rm = TRUE)
mean(boot_K80, na.rm = TRUE)
mean(boot_TN, na.rm = TRUE)
mean(boot_LOG, na.rm = TRUE)

# Conversion vers format phangorn
dna_phy <- phyDat(as.matrix(alignment_default), type="DNA")

# Test de modèles
model_test <- modelTest(dna_phy)

model_test

# Modèle avec le plus petit AICc
best_model <- model_test$Model[which.min(model_test$AICc)]
best_model

# arbre NJ initial
tree_start <- nj(dist.dna(dna))

# modèle initial
fit <- pml(tree_start, data=dna_phy)

# optimisation du modèle (exemple pour GTR+G+I)
fit_opt <- optim.pml(
    fit,
    model="GTR",
    optInv=TRUE,
    optGamma=TRUE
)

# arbre final
plot(fit_opt$tree, main=paste("ML tree -", best_model), cex=0.8)

# bootstrap ML
bs <- bootstrap.pml(fit_opt, bs=1000)

plotBS(fit_opt$tree, bs)


#################################################
# ------ d) Alignement des acides aminés ------ #
#################################################

# Choisir l'alignement à tester : cadre de lecture 1 et cadre modifié 2
aa_alignments <- list(
    frame1 = aa_alignment_1,
    frame2 = aa_alignment_2
)

# Modèles d'évolution des protéines
aa_models <- c("LG", "JTT", "Blosum62")

results <- list()

par(mfrow = c(2, 3), mar = c(2, 2, 4, 1))

# Conversion en phyDat (phangorn)
aa_phy <- lapply(aa_alignments, function(x) phyDat(as.matrix(x), type="AA"))
# Liste pour stocker les résultats

for (frame_name in names(aa_phy)) {
    phy_data <- aa_phy[[frame_name]]

    for (model in aa_models) {
        cat("\n--- Frame:", frame_name, "- Model:", model, "---\n")

        # Calcul distance et arbre
        dist_aa <- dist.ml(phy_data, model = model)
        tree_nj <- NJ(dist_aa)

        # Bootstrap (Remettre bs=1000 pour le rapport final)
        bs_trees <- bootstrap.phyDat(phy_data,
                                     FUN = function(x) NJ(dist.ml(x, model = model)),
                                     bs = 10)

        # Intégrer les valeurs de bootstrap à l'arbre original
        # p = 0 permet de garder tous les labels pour le calcul de la moyenne
        tree_with_bs <- plotBS(tree_nj, bs_trees, type = "none", p = 0)

        # Stockage
        results[[paste(frame_name, model, sep="_")]] <- tree_with_bs

        # Affichage
        plot(tree_with_bs, main=paste("NJ -", model, "-", frame_name), cex=0.8)

        # Calcul de la moyenne robuste
        # On extrait les labels, on retire les vides, on convertit en chiffre
        boot_values <- as.numeric(tree_with_bs$node.label)
        mean_bs <- mean(boot_values, na.rm = TRUE)

        cat("Robustesse moyenne bootstrap:", round(mean_bs, 2), "%\n")
    }
}

# --- RÉSUMÉ DE LA ROBUSTESSE ---
cat("\n--- RÉSUMÉ DE LA ROBUSTESSE ---\n")

# Création d'un data.frame à partir de la liste 'results'
df_summary <- data.frame(
    Analyse = names(results),
    Mean_Bootstrap = sapply(results, function(x) {
        valeurs <- as.numeric(x$node.label)
        mean(valeurs, na.rm = TRUE)
    })
)

print(df_summary)

# Trouver le gagnant pour le cadre 1 (les 3 premières lignes)
best_model_idx <- which.max(df_summary$Mean_Bootstrap[1:3])
cat("\nModèle le plus robuste (Frame 1) :", df_summary$Analyse[best_model_idx], "\n")

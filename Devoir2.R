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
seqs <- readDNAStringSet("Cetacea COI.fasta", format = "fasta")

# Alignement avec paramètres par défaut
alignment_default <- AlignSeqs(seqs)

# Alignement avec pénalités d’indels = 0
alignment_no_gap <- AlignSeqs(seqs, gapOpening = 0, gapExtension = 0)

# Comparer les alignements
width(alignment_default)
width(alignment_no_gap)

# Nombre total de gaps
sum(letterFrequency(alignment_default, "-"))
sum(letterFrequency(alignment_no_gap, "-"))

# Export pour vérification
writeXStringSet(alignment_default, "alignment_default.fasta")
writeXStringSet(alignment_no_gap, "alignment_no_gap.fasta")

##################################################
# ----- b) Traduction et cadre de lecture ------ #
##################################################

# Alignement d'acides aminées
aa_alignment_NA <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = NA)
writeXStringSet(aa_alignment_NA, "aa_alignment_NA.fasta")

aa_alignment_1 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = 1)
writeXStringSet(aa_alignment_1, "aa_alignment_1.fasta")

# Vérification
# Nombre de codons stop par séquence
stops_NA <- letterFrequency(aa_alignment_NA, "*")
stops_frame1 <- letterFrequency(aa_alignment_1, "*")

# Résumé statistique (min, médiane, max)
summary(stops_NA)
summary(stops_frame1)

# Distribution : combien de séquences ont 0, 1, 2 stops, etc.
table(stops_NA)
table(stops_frame1)

# Identifier les séquences contenant au moins un codon stop
names(aa_alignment_NA)[stops_NA > 0]
names(aa_alignment_1)[stops_frame1 > 0]


#############################################
# ------ c) Phylogénie nucléotidique ------ #
#############################################

# Conversion de l'alignement en format DNAbin (ape)
dna <- as.DNAbin(alignment_default)

# Matrices de distances avec les différents modèles
dist_JC  <- dist.dna(dna, pairwise.deletion = FALSE, model = "JC69")    # Jukes-Cantor
dist_K80 <- dist.dna(dna, pairwise.deletion = FALSE, model = "K80")     # Kimura 2 paramètres
dist_TN  <- dist.dna(dna, pairwise.deletion = FALSE, model = "TN93")    # Tamura-Nei
dist_GG95 <- dist.dna(dna, pairwise.deletion = FALSE, model = "GG95")  # Galtier-Gouy (LogDet)

# Construction des arbres
arbre_JC  <- nj(dist_JC)
arbre_K80 <- nj(dist_K80)
arbre_TN  <- nj(dist_TN)
arbre_GG95 <- nj(dist_GG95)

par(mfrow=c(2,2), mar=c(1,1,3,1))
plot(arbre_JC, cex = 0.6, main="Jukes et Cantor")
plot(arbre_K80, cex = 0.6, main="Kimura 2 parametres 1980")
plot(arbre_TN, cex = 0.6, main="Tamura et Nei 1993")
plot(arbre_GG95, cex = 0.6, main="Galtier and Gouy 1995")

#Comparaison des modèles entre eux par corrélation (Pearson)
round(cor(cbind(dist_JC, dist_TN,dist_K80, dist_GG95)),3)

# boot.phylo nécessite une matrice
dna_matrix <- as.matrix(dna)

# Bootstrap (1000 itérations)
boot_JC <- boot.phylo(phy = arbre_JC, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "JC69")), B = 1000)
boot_K80 <- boot.phylo(phy = arbre_K80, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "K80")), B = 1000)
boot_TN <- boot.phylo(phy = arbre_TN, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "TN93")), B = 1000)
boot_GG95 <- boot.phylo(phy = arbre_GG95, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "GG95")), B = 1000)

summary(boot_JC)
summary(boot_K80)
summary(boot_TN)
summary(boot_GG95)

par(mfrow=c(2,2), mar=c(1,1,3,1))
plot(arbre_JC, main="Jukes et Cantor")
nodelabels(boot_JC/10, frame="circle", bg="white", cex=0.6)

plot(arbre_K80, main="Kimura 2 paramètres")
nodelabels(boot_K80/10, frame="circle", bg="white", cex=0.6)

plot(arbre_TN, main="Tamura-Nei")
nodelabels(boot_TN/10, frame="circle", bg="white", cex=0.6)

plot(arbre_GG95, main="Galtier-Gouy")
nodelabels(boot_GG95/10, frame="circle", bg="white", cex=0.6)

round(cor(cbind(boot_JC, boot_TN, boot_K80, boot_GG95),
          use = "pairwise.complete.obs"),3)

# Conversion vers format phangorn
dna_phy <- as.phyDat(dna)

# Test de modèles
model_test <- phangorn::modelTest(dna_phy)
model_test

# tester le smodèles évolutifs
env <- attr(model_test, "env")
ls(env = env)

best_model <- model_test[which.min(model_test$AIC), ]
best_model

dm <- dist.ml(dna_phy)

treeNJ <- NJ(dm)

plot(treeNJ, main="Neighbor-Joining")

fit <- pml(treeNJ, data = dna_phy)

fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(
    fitGTR,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
)

bs <- bootstrap.pml(
    fitGTR,
    bs = 1000,
    optNni = TRUE,
    control = pml.control(trace = 0)
)

plotBS(
    midpoint(fitGTR$tree),
    bs,
    p = 0,
    type = "p",
    frame = "circle",
    cex = 0.6,
    bs.adj = c(0.5,0.5),
    bg = "white"
)






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


###########################################################
# ------ e) Analyse modelTest par cadre de lecture ------ #
###########################################################

# Création des 3 cadres de lecture (Nucléotides)
dna_frames <- list(
    frame1 = seqs,
    frame2 = subseq(seqs, start = 2),
    frame3 = subseq(seqs, start = 3)
)

results_ml_frames <- list()

# Boucle sur les 3 cadres
for (f_name in names(dna_frames)) {
    cat("\n\n=== ANALYSE :", f_name, "===\n")

    # Alignement et conversion
    aln <- AlignSeqs(dna_frames[[f_name]], verbose = FALSE)
    phy_dat <- phyDat(as.matrix(aln), type = "DNA")

    # i- Trouver le meilleur modèle avec modelTest()
    cat("Recherche du meilleur modèle...\n")
    mt <- modelTest(phy_dat)
    best_mod <- mt$Model[which.min(mt$AICc)]
    cat("Meilleur modèle (AICc) pour", f_name, ":", best_mod, "\n")

    # Reconstruction ML avec le modèle trouvé
    # On part d'un arbre NJ pour l'optimisation
    tree_init <- nj(dist.dna(as.DNAbin(phy_dat), model = "TN93"))
    fit <- pml(tree_init, data = phy_dat)

    # Mise à jour du modèle selon modelTest (on simplifie vers le modèle de base)
    # parse_model est une simplification, on utilise souvent GTR ou le modèle spécifique
    fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE)

    # Bootstrap ML (1000 itérations)
    cat("Calcul du Bootstrap ML (1000 itérations)...\n")
    bs_ml <- bootstrap.pml(fit_opt, bs = 1000, optInv = TRUE, optGamma = TRUE)

    # Sauvegarde et affichage
    results_ml_frames[[f_name]] <- list(fit = fit_opt, bs = bs_ml, model = best_mod)

    plotBS(fit_opt$tree, bs_ml, main = paste("ML Tree -", f_name, "-", best_mod))
}

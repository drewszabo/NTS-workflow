library(patRoon)
library(MS2Quant)
library(MS2Tox)
library(tidyverse)


#-----------------
# Summary Stats
#-----------------

# Number of features per sample
# Calculate the sum of non-zero values in each column of the 'resultsfGroups' data frame.
colSums(resultsfGroupsNTS != 0)

# Number of samples per feature
# Calculate the sum of non-zero values in columns 4 through 16 of the 'resultsfGroups' data frame.
samplesPerFeature <- rowSums(resultsfGroupsNTS[,4:16] != 0)

# Remove values that are equal to 0 from 'samplesPerFeature'.
samplesPerFeature <- samplesPerFeature[samplesPerFeature != 0]

# Create a histogram of the 'samplesPerFeature' values with frequency counts.
featDist <- hist(samplesPerFeature, freq = TRUE)

# Fit a linear model to the density values of the first 13 points in the histogram.
lm(featDist$density[1:13] ~ featDist$breaks)

# Create a plot of the logarithm (base 10) of the 'featDist$breaks' values (excluding the first and last points)
# against the 'featDist$density' values (excluding the first point).
plot(log10(featDist$breaks[2:12]), featDist$density[2:12])

# -------------------------
# MS2Tox - Toxicity Predictions (patRoon)
# -------------------------

compoundsSIR <- predictTox(compoundsSIR, type = "FP", concUnit = "mM")

fGroupsHaz <- patRoon::filter(fGroupsNTS, results = compoundsSIR)

fGroupsHaz <- calculateTox(fGroupsHaz, compoundsSIR)

# Save Workspace
save.image(".RData")

# -------------------------
# MS2Quant - quantification prediction
# -------------------------

eluent <- data.frame(
  time = c(0, 1200, 1500, 1560, 1800),
  B = c(5, 100, 100, 5, 5))

write.csv(eluent, "eluent.csv")

calibrants <- read.csv("calibrants.csv")

# Sometimes error in rJava when calculating rcdk::get.fingerprint. Close and remove history and project cache folder works.

compoundsSIR <- predictRespFactors(
  compoundsSIR,
  fGroupsNTS,
  calibrants = calibrants,
  eluent = eluent,
  organicModifier = "MeCN",
  pHAq = 2.7,
  concUnit = "M",
  calibConcUnit = "M",
  type = "FP")

fGroupsHaz <- calculateConcs(fGroupsHaz,
                             featureAnn = compoundsSIR,
                             areas = TRUE)

resultsfGroupsHaz <- patRoon::as.data.table(fGroupsHaz,
                                            area = TRUE,
                                            toxAggrParams = getDefPredAggrParams(),
                                            concAggrParams = getDefPredAggrParams())

resultsExport <- resultsfGroupsHaz %>%
  select(-contains("-P")) %>% # remove pooled samples
  rowwise() %>%
  mutate(maxConc = max(across(contains("_conc")))) %>%
  ungroup() %>%
  mutate(priorityScore = maxConc/(LC50/1000)) %>%
  select(group, ret, mz, adduct, maxConc, LC50, priorityScore) %>%
  left_join(resultsMB, by = "group") %>%
  left_join(resultsMoNA, by = "group") %>%
  left_join(resultsSIR, by = "group") %>%
  left_join(resultsMF, by = "group")

# Tanimoto similarity

tan.sim <- data.frame(group = character(),
                      tanimoto.sim = numeric())

for(i in 1:length(resultsExport$SMILESSIR)){
  group <- resultsExport$group[i]
  smilesSIR <- resultsExport$SMILESSIR[i]
  smilesMF <- resultsExport$SMILESMF[i]
  
  if(!is.na(smilesSIR) & !is.na(smilesMF)){
    
    sdfSIR <- smiles2sdf(smilesSIR)
    apSIR <- sdf2ap(sdfSIR)
    
    sdfMF <- smiles2sdf(smilesMF)
    apMF <- sdf2ap(sdfMF)
    
    tanimoto.sim <- cmp.similarity(apSIR, apMF)
    
  } else tanimoto.sim <- NA
  
  tan.sim <- tan.sim %>%
    add_row(group = group,
            tanimoto.sim = tanimoto.sim)
  
}


# Fileter and export results

resultsPriority <- resultsExport %>%
  left_join(tan.sim, by = "group") %>%
  dplyr::filter(scoreMB > 0 | scoreMoNA > 0 | scoreSIR > -50 | scoreMF > 5.5 | tanimoto.sim > 0.7) %>%
  mutate(silicoMatch = case_when(formulaSIR == formulaMF ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  arrange(desc(priorityScore))

write.csv(resultsPriority, "resultsPriority.csv")

# Save Workspace
save.image(".RData")

# Chromatogram of results

for(i in 1:20){
  
  png(file = paste0("chrom/priority_", i, ".png"),
      width = 400, height = 400)
  
  plotChroms(fGroups[,resultsPriority$group[i]],
             retMin = TRUE,
             showPeakArea = TRUE,
             showFGroupRect = FALSE,
             colourBy = "fGroups",
             EICParams = getDefEICParams(topMost = 1),
             showLegend = FALSE,
             xlim = c(15,25),
             title = paste0(case_when(!is.na(resultsPriority$compoundNameMB[i]) ~ resultsPriority$compoundNameMB[i],
                                      resultsPriority$compoundNameSIR[i] == "null" ~ resultsPriority$compoundNameMF[i],
                                      TRUE ~ resultsPriority$compoundNameSIR[i]),
                            "\n m/z = ",
                            round(resultsPriority$mz[i], 4),
                            ", Rt = ",
                            round(resultsPriority$ret[i]/60, 2)))
  
  dev.off()
  
}

plotChroms(fGroups[,resultsPriority$group[1]],
           retMin = TRUE,
           showPeakArea = TRUE,
           showFGroupRect = FALSE,
           colourBy = "fGroups",
           EICParams = getDefEICParams(topMost = 1),
           title = paste0(resultsPriority$compoundNameSIR[1], "\n m/z = ", round(resultsPriority$mz[1], 4), ", Rt = ", round(resultsPriority$ret[1]/60, 2)))

#-------------------------
# Level 1 confirmation (SU - KLARA)
#-------------------------

klara <- readRDS("C:/Users/drsz9242/OneDrive - Kruvelab/Drew Szabo/R/Priority Lists/Klara/klara_clean.rds")

resultsMS2Tox <- resultsMS2Tox %>%
  mutate(klara = case_when(InChIKeyMB %in% klara$StdInChIKey ~ TRUE,
                           InChIKeyMF %in% klara$StdInChIKey ~ TRUE,
                           SMILESSIR %in% klara$SMILES ~ TRUE,
                           TRUE ~ NA))

write.csv(resultsMS2Tox, "resultsMS2Tox.csv")


#-----------------
# Fold Change Analysis
#-----------------

# Using all features from all groups
myFCParams <- getFCParams(rGroups = c("HM02B", "HM03W"), # Groups for comparison
                          thresholdFC = 0.5,
                          thresholdPV = 0.05)

as.data.table(fGroupsNTS, FCParams = myFCParams)[, c("group", "FC", "FC_log", "PV", "PV_log", "classification")]

tab <- as.data.table(fGroupsNTS, FCParams = myFCParams)

# Volcano plot
plotVolcano(fGroupsNTS, myFCParams,
            xlim = c(-15, 15))

# only keep feature groups that significantly increase or decrease
fGroupsFC <- fGroupsNTS[, tab[classification %in% c("increase", "decrease")]$group]

# Prepare report
resultsfGroupsFC <- patRoon::as.data.table(fGroupsFC, area = TRUE, average = TRUE)

resultsfGroupsFC <- resultsfGroupsFC %>%
  left_join(resultsMS2Tox, by = "group") %>%
  left_join(MB_SS, by = "group") %>%
  left_join(SIR_SS, by = "group") %>%
  left_join(MF_SS, by = "group") %>%
  dplyr::filter(HM02B != 0 & HM03W != 0,
                scoreMB > 0.8 | scoreSIR > -100| scoreMF > 3)



#-----------------
# PCA Analysis
#-----------------

library(tidyverse)
library(factoextra)

# Select only the sample abundances
da <- resultsfGroups %>%
  column_to_rownames("group") %>%
  select(3:15)

# Remove rows with all zero values (from QC removal)
da <- da[rowSums(da[])>0,]

# Remove sample identifiers
colnames(da) <- gsub("HM", "", colnames(da))

# PCA
results.pca <- prcomp(da, scale. = TRUE)

# Scree plot
fviz_eig(results.pca)

# Feature plot
fviz_pca_ind(results.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Sample plot
fviz_pca_var(results.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Biplot
biplot.pca <- fviz_pca_biplot(
  results.pca,
  label = "var",
  col.ind = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
)
biplot.pca

factors.pca <- data.frame(group = biplot.pca$data$name,
                          pc1 = biplot.pca$data$x,
                          pc2 = biplot.pca$data$y,
                          cos2 = biplot.pca$data$cos2)

factors.pca <- factors.pca %>%
  left_join(MB_SS, by = "group") %>%
  left_join(SIR_SS, by = "group") %>%
  left_join(MF_SS, by = "group") %>%
  dplyr::filter(!is.na(scoreMB) | !is.na(scoreSIR) | !is.na(scoreMB))

ggsave("nts_pca.png",
       width = 6,
       height = 6)


#-----------------
# Bray-Curtis Dendogram
#-----------------

library(tidyverse)
library(vegan)

# Bray Curtis distance
bcDist <- vegdist(da, method = "bray")
bcHclust <- hclust(bcDist)

fviz_dend(bcHclust, k = 4)

# Calculate cophenetic correlation coefficient
bcCoph <- cophenetic(bcHclust)
cor(bcDist, bcCoph)

# Euclidian distance
euclDist <- vegdist(da, method = "euclidean")
euclHclust <- hclust(euclDist)

fviz_dend(euclHclust, k = 4)

# Calculate cophenetic correlation coefficient
euclCoph <- cophenetic(euclHclust)
cor(euclDist, euclCoph)

# Chisq distance
chiDist <- vegdist(da, method = "chisq")
chiHclust <- hclust(chiDist)

fviz_dend(euclHclust, k = 6)

# Calculate cophenetic correlation coefficient
chiCoph <- cophenetic(chiHclust)
cor(chiDist, chiCoph)


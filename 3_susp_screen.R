# -------------------------
# Athens RTI prediction model
# -------------------------

# RTI Calibration
calList <- read.csv("RTI_calibrant.csv")
calList <- calList %>%
  select(name, SMILES, rt, RTI)

fGroupsCal <- screenSuspects(fGroups, calList, adduct = "[M+H]+", mzWindow = 0.005, rtWindow = 30, onlyHits = TRUE)
resultsfGroupsCal <- patRoon::as.data.table(fGroupsCal, area = TRUE, average = TRUE)

resultsfGroupsCal <- resultsfGroupsCal %>%
  rename(name = susp_name) %>%
  left_join(calList, by = "name")

ggplot(data = resultsfGroupsCal) +
  geom_point(aes(x = ret, y = RTI)) +
  geom_smooth(aes(x = ret, y = RTI),
              method = "lm",
              linetype = 2,
              colour = "black") +
  geom_text_repel(aes(x = ret, y = RTI, label = name)) +
  theme_classic()

ggsave("rti_model.png",
       width = 5,
       height = 5)

calModel <- lm(data = resultsfGroupsCal, formula = RTI ~ ret)
summary(calModel)
rtCoeff <- calModel$coefficients[[2]]
rtInt <- calModel$coefficients[[1]]

# Calculated RMSE of model
rmse <- sqrt(mean((resultsfGroupsCal$ret - calModel$fitted.values)^2))
rmse_rt <- (rmse - rtInt) / rtCoeff

#---------------------
# KEMI suspect list (textile substances)
#---------------------

kemi <- read.csv("kemiFischer_clean.csv")

kemi[kemi == ""] <- NA

kemi <- kemi %>%
  select(Name, StdInChIKey, SMILES, Molecular_Formula, Monoiso_Mass, exposureScore, Pred_RTI_Positive_ESI, Pred_RTI_Negative_ESI) %>%
  distinct(StdInChIKey, .keep_all = TRUE) %>%
  rename(name = StdInChIKey)

fGroupsKemi <- screenSuspects(fGroupsNTS, kemi, adduct = "[M+H]+", mzWindow = 0.002, onlyHits = TRUE)

# Visualise EIC
checkFeatures(fGroupsKemi, clearSession = TRUE)
fGroupsKemi <- patRoon::filter(fGroupsKemi, checkFeaturesSession = "checked-features.yml")

# Prepare SS Report
resultsfGroupsKemi <- patRoon::as.data.table(fGroupsKemi, area = TRUE, average = TRUE)

# Number of hits per features
susp_hits <- str_count(resultsfGroupsKemi$susp_name, ",") + 1
sum(susp_hits)
min(susp_hits)
max(susp_hits)
mean(susp_hits)
susp_dist <- hist(susp_hits, freq = FALSE)
sum(susp_dist$counts[2:length(susp_dist$counts)])

# Pivot multiple susp hits to rows
resultsfGroupsKemi <- resultsfGroupsKemi %>%
  separate_wider_delim(susp_name, delim = ",", names = as.character(seq(1, max(susp_hits))), too_few = "align_start", too_many = "error") %>%
  gather(key = "susp_number", value = "name", -colnames(resultsfGroupsKemi)[1:18]) %>%
  drop_na(name)

# Match susp with table (RTI, name, etc)
resultsfGroupsKemi <- resultsfGroupsKemi %>%
  left_join(kemi, by = "name") %>%
  rename(StdInChIKey = name)

# Apply and filter RTI
resultsfGroupsKemi <- resultsfGroupsKemi %>%
  mutate(pred_ret = (Pred_RTI_Positive_ESI - rtInt) / rtCoeff,
         rt_filter = ifelse(pred_ret - ret < rmse_rt & pred_ret - ret > -rmse_rt, TRUE, FALSE))

resultsKemi <- resultsfGroupsKemi %>%
  left_join(resultsMB, by = "group") %>%
  left_join(resultsSIR, by = "group") %>%
  left_join(resultsMF, by = "group") %>%
  left_join(resultsMoNA, by = "group") %>%
  dplyr::filter(rt_filter %in% c(NA, TRUE))

# Tanimoto similarity

tan.sim <- data.frame(group = character(),
                      StdInChIKey = character(),
                      matchMB = numeric(),
                      matchSIR = numeric(),
                      matchMF = numeric())

for(i in 1:length(resultsKemi$StdInChIKey)){
  group <- resultsKemi$group[i]
  StdInChIKey <- resultsKemi$StdInChIKey[i]
  
  # MassBank
  if(!is.na(resultsKemi$SMILES[i]) & !is.na(resultsKemi$SMILESMB[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsKemi$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfMB <- ChemmineR::smiles2sdf(resultsKemi$SMILESMB[i])
    apMB <- ChemmineR::sdf2ap(sdfMB)
    
    matchMB <- ChemmineR::cmp.similarity(apSS, apMB)
    
  } else matchMB <- NA
  
  # SIRIUS
  if(!is.na(resultsKemi$SMILES[i]) & !is.na(resultsKemi$SMILESSIR[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsKemi$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfSIR <- ChemmineR::smiles2sdf(resultsKemi$SMILESSIR[i])
    apSIR <- ChemmineR::sdf2ap(sdfSIR)
    
    matchSIR <- ChemmineR::cmp.similarity(apSS, apSIR)
    
  } else matchSIR <- NA
  
  # MetFrag
  if(!is.na(resultsKemi$SMILES[i]) & !is.na(resultsKemi$SMILESMF[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsKemi$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfMF <- ChemmineR::smiles2sdf(resultsKemi$SMILESMF[i])
    apMF <- ChemmineR::sdf2ap(sdfMF)
    
    matchMF <- ChemmineR::cmp.similarity(apSS, apMF)
    
  } else matchMF <- NA
  
  tan.sim <- tan.sim %>%
    add_row(group = group,
            StdInChIKey = StdInChIKey,
            matchMB = matchMB,
            matchSIR = matchSIR,
            matchMF = matchMF)
  
}

resultsKemi <- resultsKemi %>%
  left_join(tan.sim, by = c("group", "StdInChIKey"))

write.csv(resultsKemi, "resultsKemiID_pos.csv")



# Plot chromatograms matched with structural annotation
plotChroms(fGroupsKemi[, resultsSS$group],
           EICParams = getDefEICParams(topMost = 1), # only most intense feature in each group
           showPeakArea = TRUE, # show integrated areas
           showFGroupRect = FALSE,
           showLegend = FALSE,
           retMin = TRUE,
           title = "")


ggplot(data = resultsfGroupsKemi) +
  geom_point(aes(x = ret/60,
                 y = pred_ret/60),
             alpha = ifelse(resultsfGroupsKemi$rt_filter == TRUE, 1, 0.1)) +
  geom_abline(slope = 1) +
  scale_x_continuous(name = "Measured RT (min)") +
  scale_y_continuous(name = "Predicted RT (min)") +
  theme_classic()

ggsave("susp_screen_RTI.png",
       width = 5,
       height = 5)

# Save Workspace
save.image(".RData")

#-----------------------
# REACH suspect list
#-----------------------

reach <- read.csv("REACH_PMT.csv")

reach[reach == ""] <- NA

reach <- reach %>%
  drop_na(StdInChIKey) %>%
  distinct(StdInChIKey, .keep_all = TRUE) %>%
  rename(name = StdInChIKey,
         formula = Molecular_Formula)
  
reach.ss <- reach %>%  
  select(name, formula)

fGroupsReach <- screenSuspects(fGroupsNTS, reach.ss, adduct = "[M+H]+", mzWindow = 0.002, onlyHits = TRUE)

# Visualise EIC
checkFeatures(fGroupsReach, clearSession = TRUE)
fGroupsReach <- patRoon::filter(fGroupsReach, checkFeaturesSession = "checked-features.yml")

# Prepare SS Report
resultsfGroupsReach <- patRoon::as.data.table(fGroupsReach, area = TRUE, average = TRUE)

# Number of hits per features
susp_hits <- str_count(resultsfGroupsReach$susp_name, ",") + 1
min(susp_hits)
max(susp_hits)
mean(susp_hits)
susp_dist <- hist(susp_hits, freq = FALSE)
sum(susp_dist$counts[2:length(susp_dist$counts)])

# Pivot multiple susp hits to rows
resultsfGroupsReach <- resultsfGroupsReach %>%
  separate_wider_delim(susp_name, delim = ",", names = as.character(seq(1, max(susp_hits))), too_few = "align_start", too_many = "error") %>%
  gather(key = "susp_number", value = "name", -colnames(resultsfGroupsReach)[1:18]) %>%
  drop_na(name)

# Match susp with table (RTI, name, etc)
resultsfGroupsReach <- resultsfGroupsReach %>%
  left_join(reach, by = "name") %>%
  rename(StdInChIKey = name)

# Apply and filter RTI
resultsfGroupsReach <- resultsfGroupsReach %>%
  mutate(pred_ret = (Pred_RTI_Positive_ESI - rtInt) / rtCoeff,
         rt_filter = ifelse(pred_ret - ret < rmse_rt & pred_ret - ret > -rmse_rt, TRUE, FALSE))

resultsReach <- resultsfGroupsReach %>%
  left_join(resultsMB, by = "group") %>%
  left_join(resultsSIR, by = "group") %>%
  left_join(resultsMF, by = "group") %>%
  left_join(resultsMoNA, by = "group") %>%
  dplyr::filter(rt_filter %in% c(NA, TRUE))

# Tanimoto similarity

tan.sim <- data.frame(StdInChIKey = character(),
                      matchMB = numeric(),
                      matchSIR = numeric(),
                      matchMF = numeric())

for(i in 1:length(resultsReach$StdInChIKey)){
  StdInChIKey <- resultsReach$StdInChIKey[i]
  
  # MassBank
  if(!is.na(resultsReach$SMILES[i]) & !is.na(resultsReach$SMILESMB[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsReach$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfMB <- ChemmineR::smiles2sdf(resultsReach$SMILESMB[i])
    apMB <- ChemmineR::sdf2ap(sdfMB)
    
    matchMB <- ChemmineR::cmp.similarity(apSS, apMB)
    
  } else matchMB <- NA
  
  # SIRIUS
  if(!is.na(resultsReach$SMILES[i]) & !is.na(resultsReach$SMILESSIR[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsReach$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfSIR <- ChemmineR::smiles2sdf(resultsReach$SMILESSIR[i])
    apSIR <- ChemmineR::sdf2ap(sdfSIR)
    
    matchSIR <- ChemmineR::cmp.similarity(apSS, apSIR)
    
  } else matchSIR <- NA
  
  # MetFrag
  if(!is.na(resultsReach$SMILES[i]) & !is.na(resultsReach$SMILESMF[i])){
    
    sdfSS <- ChemmineR::smiles2sdf(resultsReach$SMILES[i])
    apSS <- ChemmineR::sdf2ap(sdfSS)
    
    sdfMF <- ChemmineR::smiles2sdf(resultsReach$SMILESMF[i])
    apMF <- ChemmineR::sdf2ap(sdfMF)
    
    matchMF <- ChemmineR::cmp.similarity(apSS, apMF)
    
  } else matchMF <- NA
  
  tan.sim <- tan.sim %>%
    add_row(StdInChIKey = StdInChIKey,
            matchMB = matchMB,
            matchSIR = matchSIR,
            matchMF = matchMF)
  
}

resultsReach <- resultsReach %>%
  left_join(tan.sim, by = "StdInChIKey")

write.csv(resultsReach, "resultsReachID_pos.csv")


# Plot chromatograms matched with structural annotation

for(i in 1:length(resultsReach$group)){
  
  png(filename = paste0("reach/", resultsReach$group[i], ".png"))
  
  plotChroms(fGroupsReach[, resultsReach$group[i]],
             colourBy = "fGroups",
             EICParams = getDefEICParams(topMost = 1), # only most intense feature in each group
             showPeakArea = TRUE, # show integrated areas
             showFGroupRect = FALSE,
             showLegend = FALSE,
             retMin = TRUE,
             title = "",
             xlim = c((resultsReach$ret.x[i] - 180)/60, (resultsReach$ret.x[i] + 180)/60))
  
  dev.off()
  
}

plotChroms(fGroupsReach[, resultsReach$group],
           colourBy = "fGroups",
           EICParams = getDefEICParams(topMost = 1), # only most intense feature in each group
           showPeakArea = TRUE, # show integrated areas
           showFGroupRect = FALSE,
           showLegend = FALSE,
           retMin = TRUE,
           title = "")

# Plot chromatograms matched with structural annotation
plotChroms(fGroupsReach[, resultsSSReach$group[2]],
           colourBy = "rGroup",
           showPeakArea = TRUE, # show integrated areas
           showFGroupRect = FALSE,
           showLegend = FALSE,
           retMin = TRUE,
           title = "")


ggplot(data = resultsfGroupsReach) +
  geom_point(aes(x = ret/60,
                 y = pred_ret/60),
             alpha = ifelse(resultsfGroupsReach$rt_filter == TRUE, 1, 0.1)) +
  geom_abline(slope = 1) +
  scale_x_continuous(name = "Measured RT (min)") +
  scale_y_continuous(name = "Predicted RT (min)") +
  theme_classic()

ggsave("susp_screen_RTI.png",
       width = 5,
       height = 5)

# Save Workspace
save.image(".RData")


#-----------------------------
# ClassyFire annotation
#-----------------------------

# Only works for InChIKey for now
getClassy <- function(screenTable){
  
  for (i in 1:length(screenTable$StdInChIKey)) {
    
    # Get classification of compound
    classification <- get_classification(screenTable$StdInChIKey[i])
    nLevels <- length(classification(classification)[[1]]) # number of levels
    
    # class level (almost always available)
    if(nLevels > 2) {
    className <- classification(classification)[[3,2]] # level 3
    } else className <- NA
    
    # subclass level if available
    if(nLevels > 3) {
      subclassName <- classification(classification)[[4,2]] # level 4
      
    } else subclassName <- NA
    
    screenTable$class[i] <- className
    screenTable$subclass[i] <- subclassName
    
  }
  
  return(screenTable)
  
}

resultsSSReach <- getClassy(resultsSSReach)
resultsSS <- getClassy(resultsSS)
resultsMB <- getClassy(resultsMB %>% rename(StdInChIKey = InChIKey))



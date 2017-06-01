# LOAD IN FILE

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/CORE_Forbotanycluster/COREOTUs_ExNWater")
#################### LOAD FILES #####################
coreData <- read.delim("CORE_OTUS.txt"
                       , stringsAsFactors = FALSE
                       , header = FALSE
                       , row.names = 1)
abundanceData <- read.delim("OTU_Table_relabund.txt"
                            , stringsAsFactors = FALSE
                            , header = TRUE
                            , row.names = 1)
metadata <- read.delim("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/LFSE2/WHOLE_OTU_OTHER/MF_for_Lfse.txt"
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1)
taxaID <- read.delim("taxaIDLegend.txt"
                     , stringsAsFactors = FALSE
                     , header = TRUE
                     , row.names = 1)

nullTreatment <- 'W'
treatment <- 'ExNWater'
outputGraphName <- paste0("OTUEnrichment_",treatment,".jpeg")

#################### FUNCTIONS #####################
getUniqueCore <- function(coreData) {
  allCoreOTUs <- c()
  for (treatmentType in rownames(coreData)) {
    allCoreOTUs <- c(allCoreOTUs, coreData[treatmentType,])
  }
  uniqueAllCoreOTUs <- unique(unlist(allCoreOTUs))
  uniqueAllCoreOTUs <- uniqueAllCoreOTUs[!is.na(uniqueAllCoreOTUs)]
  return(uniqueAllCoreOTUs)
}

getTreatmentTypes <- function(coreData = coreData, treatment = treatment, nullTreatment = nullTreatment) {
  treatmentTypes <- rownames(coreData)[grep(treatment, rownames(coreData))]
  treatmentTypes <- gsub(paste0(treatment,":"), "", treatmentTypes)
  treatmentTypes <- treatmentTypes[-grep(paste0("^",nullTreatment,"$"), treatmentTypes)]
  return(treatmentTypes)
}

getTreatmentNames <- function(treatmentTypes = treatmentTypes, metadata = metadata, treatment = treatment) {
  treatmentNames <- list()
  for (i in treatmentTypes) { 
    treatmentNames[[paste(i)]] <- c(rownames(metadata)[grep(paste0("^",i,"$"), metadata[,treatment])])
  }
  return(treatmentNames)
}

getFoldChange <- function(treatmentTypes = treatmentTypes,treatment = treatment, coreData = coreData, abundanceData = abundanceData) {
  # Find abundances of each; compare with NULL
  averageAbundances <- list()
  foldChange <- list()
  APPEARED <- list()
  for (t in treatmentTypes) {
    averageAbundances[[t]] <- list()
    foldChange[[t]] <- list()
    listOTUs <- coreData[grep(paste0(treatment,":",t), rownames(coreData)),]
    listOTUs <- listOTUs[!is.na(listOTUs)]
    APPEARED[[t]] <- c()
    for (OTU in listOTUs) {
      abundanceVector <- c()
      for (sample in treatmentNames[[t]]) {
        x <- grep(paste0("^",OTU,"$"), rownames(abundanceData))
        y <- grep(paste0("^",sample,"$"), colnames(abundanceData))
        abundanceVector <- c(abundanceVector, abundanceData[x,y])
      }
      averageAbundances[[t]][[as.character(OTU)]] <- mean(abundanceVector)
      abundanceVector <- c()
      sample <- nullTreatmentNames[1]
      for (sample in nullTreatmentNames) {
        x <- grep(paste0("^",OTU,"$"), rownames(abundanceData))
        y <- grep(paste0("^",sample,"$"), colnames(abundanceData))
        abundanceVector <- c(abundanceVector, abundanceData[x,y])
      }
      if (mean(abundanceVector) == 0) {
        APPEARED[[t]] <- c(APPEARED[[t]], OTU)
      } else if (averageAbundances[[t]][[as.character(OTU)]]/mean(abundanceVector) == 0) { 
        foldChange[[t]][[as.character(OTU)]] <- 0
        } else {
        foldChange[[t]][[as.character(OTU)]] <- log(averageAbundances[[t]][[as.character(OTU)]]/mean(abundanceVector), base = 2)
      }
    }
  }
  for (i in 1:length(foldChange)) {
    toDelete <- c()
    for (j in 1:length(foldChange[[i]])) {
      if (abs(foldChange[[i]][[j]]) < 2) {
        toDelete <- c(toDelete, j)
      }
    }
    foldChange[[i]] <- foldChange[[i]][-toDelete]
  }
  return(list(foldChange, APPEARED))
}

getRepresentation <- function(abundanceData = abundanceData, treatment = treatment, coreData = coreData, testTreatment = testTreatment, uniqueAllCoreOTUs = uniqueAllCoreOTUs) {
  outputMatrix <- matrix(ncol = length(uniqueAllCoreOTUs), nrow = 1)
  colnames(outputMatrix) <- uniqueAllCoreOTUs
  # Get OTUs and remove NAs
  # OTUsToGetRep <- coreData[grep(paste0(treatment,":",testTreatment), rownames(coreData)),]
  # OTUsToGetRep <- OTUsToGetRep[!is.na(OTUsToGetRep)]
  # Find abundance data in abundanceData
  singleLineAbund <- c()
  for (i in uniqueAllCoreOTUs) {
    singleLineAbund <- abundanceData[grep(paste0("^",i,"$"), rownames(abundanceData)),]
    treatmentNamesTemp <- getTreatmentNames(treatmentTypes = testTreatment, metadata = metadata, treatment = treatment)
    singleLineAbundTreat <- c()
    for (j in unlist(treatmentNamesTemp)) {
      singleLineAbundTreat <- c(singleLineAbundTreat, singleLineAbund[,grep(paste0("^",j,'$'), colnames(singleLineAbund))])
    }
    outputMatrix[1,paste0(i)] <- mean(singleLineAbundTreat)
  }
  return(outputMatrix)
}


#################### START #####################

# Make list of all OTUs in dataset
uniqueAllCoreOTUs <- getUniqueCore(coreData)

# Go through each OTU in each core list and see how it is enriched compared to basal
nullTreatmentNames <- rownames(metadata)[grep(paste0("^",nullTreatment,"$"), metadata[,treatment])]

#################### FOR EXNEXN #####################

# Find treatment Types
treatmentTypes <- getTreatmentTypes(coreData = coreData,treatment = treatment, nullTreatment = nullTreatment)

# Find treatment names
treatmentNames <- getTreatmentNames(treatmentTypes = treatmentTypes, metadata = metadata, treatment = treatment)

# Calculate fold change over all OTUs for each group
foldChangeOUT <- getFoldChange(treatmentTypes = treatmentTypes
                            , treatment = treatment
                            , coreData = coreData
                            , abundanceData = abundanceData)
foldChange <- foldChangeOUT[[1]]
APPEARED <- foldChangeOUT[[2]]

forBarPlot <- matrix(nrow = length(uniqueAllCoreOTUs), ncol = length(treatmentTypes))
rownames(forBarPlot) <- uniqueAllCoreOTUs
colnames(forBarPlot) <- treatmentTypes

for (i in 1:nrow(forBarPlot)) {
  for (j in 1:ncol(forBarPlot)) {
    tempValue <- foldChange[[paste0(colnames(forBarPlot)[j])]][[paste0(rownames(forBarPlot)[i])]]
    if (is.null(tempValue)) {
      next 
    } else {
      forBarPlot[i,j] <- tempValue
    }
  }
}

#transpose
forBarPlot <- t(forBarPlot)

# Get rid of things that are always NA
allNA <- c()
for (column in 1:ncol(forBarPlot)) {
  if (all(is.na(forBarPlot[,column]))) {
    allNA <- c(allNA, column)
  }
}

forBarPlot <- forBarPlot[,-allNA]


# Reorder to make it greatest to smallest
averagesAcrossOTUs <- colSums(forBarPlot, na.rm = TRUE)
forBarPlot <- forBarPlot[,order(averagesAcrossOTUs, decreasing = TRUE)]
# forBarPlotReduce <- forBarPlot[,order(forBarPlot[nrow(forBarPlot),], decreasing = TRUE)]


# # Note things that appeared and disappeared
# appearedstuff <- forBarPlot
# for (i in colnames(forBarPlot)) {
#   for (type in rownames(forBarPlot)) {
#     if (i %in% APPEARED[[type]]) {
#       appearedstuff[type,i] <- 1
#     } else if (!is.na(forBarPlot[type,i])) {
#       appearedstuff[type,i] <- NA
#     } else {
#       appearedstuff[type,i] <- NA
#     }
#   }
# }

# Make barplot that shows representation in control treatment

# abundRepInControl <- getRepresentation(abundanceData = abundanceData, treatment = treatment, coreData = coreData, testTreatment = "W", uniqueAllCoreOTUs = uniqueAllCoreOTUs)
# abundRepInControl <- abundRepInControl[,order(averagesAcrossOTUs, decreasing = TRUE)]

# PLOT 3 GRAPHS ON TOP OF EACH OTHER
# Get rid of NA's
allValuesWithoutNA <- forBarPlot[!is.na(forBarPlot)]
minValue <- min(allValuesWithoutNA)
maxValue <- max(allValuesWithoutNA)

########## GRAPHING ##############


# FOR EXNEXN; 1 control, 3 treatments
jpeg(paste0(outputGraphName), width = 1000, height = 1000)
par(fig = c(0,1,0.64,.95), mar = c(1,8,1,1))
barplot(forBarPlot[3,]
        , col = c("green")
        , border = NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Nereocystis"
)
par(fig = c(0,1,0.33,0.64), mar = c(1,8,1,1), new = TRUE)
barplot(forBarPlot[2,]
        , col = c("red")
        , border = NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Mastocarpus"
)
par(fig = c(0,1,0.02,0.33), mar = c(1,8,1,1), new = TRUE)
barplot(forBarPlot[1,]
        , col = c("brown")
        , border= NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Nereo and Mast"
)
# par(fig = c(0,1,0,0.25), mar = c(4,8,1,1), new = TRUE)
# barplot(abundRepInControl
#         , xaxt = 'n'
#         , col = 'grey'
#         , ylab = "%Abund. of OTU in Control"
#         , xlab = '')
par(fig = c(0,1,0,1), new = TRUE)
title(main = "OTU Enrichment of treatments compared to control", cex.main = 2, line = -1)
title(xlab = "OTUs", cex.lab = 2, line = 0)
par(fig = c(0,1,0,1), new = TRUE)
title(ylab = "Amount-Fold Enrichment", cex.lab = 2, line = 5)
dev.off()


# FOR EXNWater; 1 control 4 treatments
jpeg(paste0(outputGraphName), width = 1000, height = 1000)
par(fig = c(0,1,0.72,0.95), mar = c(1,8,1,1))
barplot(forBarPlot[4,]
        , col = c("brown")
        , border = NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Nereo and Mast"
)
par(fig = c(0,1,0.49,.72), mar = c(1,8,1,1), new = TRUE)
barplot(forBarPlot[3,]
        , col = c("red")
        , border = NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Mastocarpus"
)
par(fig = c(0,1,0.26,0.49), mar = c(1,8,1,1), new = TRUE)
barplot(forBarPlot[2,]
        , col = c("green")
        , border = NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Nereocystis"
)
par(fig = c(0,1,0.03,0.26), mar = c(1,8,1,1), new = TRUE)
barplot(forBarPlot[1,]
        , col = c("yellow")
        , border= NA
        , xaxt = 'n'
        , ylim = c(minValue, maxValue)
        , xlab = ''
        , ylab = "Meristem Nereocystis"
)
# par(fig = c(0,1,0,0.19), mar = c(4,8,1,1), new = TRUE)
# barplot(abundRepInControl
#         , xaxt = 'n'
#         , col = 'grey'
#         , ylab = "% Abund. of OTU in Control"
#         , xlab = '')
par(fig = c(0,1,0,1), new = TRUE)
title(main = "OTU Enrichment of treatments compared to control", cex.main = 2, line = -1)
title(xlab = "OTUs", cex.lab = 2, line = 2)
par(fig = c(0,1,0,1), new = TRUE)
title(ylab = "Amount-Fold Enrichment", cex.lab = 2, line = 5)
dev.off()

# Curious to see how many in each?
TOTALCOUNTS <- matrix(nrow = nrow(forBarPlot))
rownames(TOTALCOUNTS) <- rownames(forBarPlot)

for (i in 1:nrow(forBarPlot)) {
  TOTALCOUNTS[rownames(forBarPlot)[i],1] <- length(forBarPlot[i,][!is.na(forBarPlot[i,])])
}

# Quick test to see what things are absent from NereoMast that are present in Nereo
# ExNExN
discrepancies <- c()
for (i in 1:ncol(forBarPlot)) {
  storeTF <- is.na(forBarPlot['B',i])
  if (!storeTF) {
    if (forBarPlot['B',i] < 0) {
      storeTF <- TRUE
    } 
  }
  if (((forBarPlot['N',i] > 0) && !is.na(forBarPlot['N',i])) && (storeTF)){
    discrepancies <- c(discrepancies, i)
  }
}




# ExNWater
discrepancies <- c()
for (i in 1:ncol(forBarPlot)) {
  storeTF <- is.na(forBarPlot['NereoMast',i])
  if (!storeTF) {
    if (forBarPlot['NereoMast',i] < 0) {
      storeTF <- TRUE
    } 
  }
  if (((forBarPlot['Nereo',i] > 0) && !is.na(forBarPlot['Nereo',i])) && (storeTF)){
    discrepancies <- c(discrepancies, i)
  }
}




# FOR BOTH
discrepBarPlot <- forBarPlot[,discrepancies]
discrepBarPlot[which(is.na(discrepBarPlot))] <- 0

# Doesn't Work
# quartz()
# barplot(discrepBarPlot[2:4,]
#         , beside = TRUE
#         , col = c('green','red','brown')
#         , space = c(rep(c(rep(0,3), 2),ncol(discrepBarPlot)-1), rep(0,3))
#         )
# ncol(discrepBarPlot)

############## GET IDS FOR IMPORTANT GROUPS ##################

# FIRST, find IDs that are Enriched, Reduced, across multiple combos.

limitedEnriched <- list()
for (i in treatmentTypes) {
  limitedEnriched[[paste0(i)]] <- c()
}
limitedReduced <- list()
for (i in treatmentTypes) {
  limitedReduced[[paste0(i)]] <- c()
}
EnrichedAll <- c()
ReducedAll <- c()

for (taxa in colnames(forBarPlot)) {
  testSame <- c()
  for (treatment in treatmentTypes) {
    if (is.na(forBarPlot[treatment,taxa])) {
      testSame <- c(testSame, FALSE)
    } else if (forBarPlot[treatment,taxa] > 0) {
      testSame <- c(testSame, TRUE)
    } else {
      testSame <- c(testSame, FALSE)
    }
  }
  if (all(testSame)) {
    EnrichedAll <- c(EnrichedAll, taxa)
  } else {
    for (i in 1:length(treatmentTypes)) {
      if ((testSame[i]) && (!testSame[-i])) {
        limitedEnriched[[treatmentTypes[i]]] <- c(limitedEnriched[[treatmentTypes[i]]], taxa)
      }
    }
  }
  
  
  testDiff <- c()
  for (treatment in treatmentTypes) {
    if (is.na(forBarPlot[treatment,taxa])) {
      testDiff <- c(testDiff, FALSE)
    } else if (forBarPlot[treatment,taxa] < 0) {
      testDiff <- c(testDiff, TRUE)
    } else {
      testDiff <- c(testDiff, FALSE)
    }
  }
  if (all(testDiff)) {
    ReducedAll <- c(ReducedAll, taxa)
  } else {
    for (i in 1:length(treatmentTypes)) {
      if ((testDiff[i]) && (!testDiff[-i])) {
        limitedReduced[[treatmentTypes[i]]] <- c(limitedReduced[[treatmentTypes[i]]], taxa)
      }
    }
  }
}

# ENRICHED

EnrichedAllTaxa <- c()
for (i in EnrichedAll) {
  EnrichedAllTaxa <- c(EnrichedAllTaxa, taxaID[grep(paste0("^", i, "$"), rownames(taxaID)),])
}

# REDUCED

ReducedAllTaxa <- c()
for (i in ReducedAll) {
  ReducedAllTaxa <- c(ReducedAllTaxa, taxaID[grep(paste0("^", i, "$"), rownames(taxaID)),])
}

# LIMITED ENRICHED

limitedEnrichedTaxa <- limitedEnriched
for (i in treatmentTypes) {
  for (j in 1:length(limitedEnriched[[i]])) {
    limitedEnrichedTaxa[[i]][j] <- taxaID[grep(paste0("^", limitedEnriched[[i]][j], "$"), rownames(taxaID)),]
  }
}
 
# LIMITED REDUCED

limitedReducedTaxa <- limitedReduced
for (i in treatmentTypes) {
  for (j in 1:length(limitedReduced[[i]])) {
    limitedReducedTaxa[[i]][j] <- taxaID[grep(paste0("^", limitedReduced[[i]][j], "$"), rownames(taxaID)),]
  }
}



####### PLOT ENRICHMENT/REDUCTION IN TAXASUMMARIES #######
ALLSamples <- as.vector(unlist(c(treatmentNames, nullTreatmentNames)))
# Make an OTU table with only desired samples
filteredAbundanceData <- matrix(nrow = nrow(abundanceData), ncol = length(unlist(ALLSamples)))
rownames(filteredAbundanceData) <- rownames(abundanceData)
colnames(filteredAbundanceData) <- unlist(ALLSamples)

for (j in 1:length(ALLSamples)) {
  # print(head(abundanceData[,grep(j, colnames(abundanceData))]))
  filteredAbundanceData[,j] <- abundanceData[,grep(ALLSamples[j], colnames(abundanceData))]
}

# Get rid of zero's
filteredAbundanceData <- filteredAbundanceData[-which(rowSums(filteredAbundanceData) == 0),]


# Sample colours Randomly
coloursRandom <- sample(colors()[grep('gr(e|a)y|white', colors(), invert = TRUE)], length(uniqueAllCoreOTUs), replace = TRUE)

colourLegendEnriched <- matrix(nrow = nrow(filteredAbundanceData))
colourLegendReduced <- matrix(nrow = nrow(filteredAbundanceData))

rownames(colourLegendEnriched) <- rownames(filteredAbundanceData)
rownames(colourLegendReduced) <- rownames(filteredAbundanceData)
count = 0
for (i in 1:length(rownames(filteredAbundanceData))) {
  taxaTemp <- rownames(filteredAbundanceData)[i]
  if (taxaTemp %in% EnrichedAll) {
    count = count + 1
    colourLegendEnriched[i,] <- coloursRandom[count]
  } else {
    colourLegendEnriched[i,] <- NA
    } 
  if (taxaTemp %in% ReducedAll) {
    count = count + 1
    colourLegendReduced[i,] <- coloursRandom[count]
  } else {
    colourLegendReduced[i,] <- NA
  }
}



# Sample labels-- will have to change each time, depending on samples
#ExNExN
sampleLab <- c(rep("Nereo+Mast", 5), rep("Mastocarpus", 4), rep("Nereocystis", 5), rep("None", 4))

# Get colour labels
OTUIDsENRICH <- rownames(colourLegendEnriched)[!is.na(colourLegendEnriched)]
OTUTAXONENRICH <- sapply(OTUIDsENRICH, function(x) {
  taxaID[grep(paste0("^",x,'$'), rownames(taxaID)),]
})
coloursLegendENRICH <- colourLegendEnriched[!is.na(colourLegendEnriched)]

OTUIDsREDUCE <- rownames(colourLegendReduced)[!is.na(colourLegendReduced)]
OTUTAXONREDUCE <- sapply(OTUIDsREDUCE, function(x) {
  taxaID[grep(paste0("^",x,'$'), rownames(taxaID)),]
})
coloursLegendREDUCE <- colourLegendReduced[!is.na(colourLegendReduced)]


jpeg("TaxaSummaries_ExNExN_Enriched.jpeg", width = 2000, height = 2000)
par(fig = c(0,1,0.5,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = colourLegendEnriched
        , las = 2
        , border = NA
        , names.arg = sampleLab
        , space = c(0,0,0,0,0,1,0,0,0,1,0,0,0,0,2,0,0,0)
        , ylab = "Relative Abundance"
        )
title(xlab = "Treatment", line = 6)
par(fig = c(0,1,0,0.5), mar = c(0,0,0,0), new = TRUE)
legend("center"
       , legend = OTUTAXONENRICH
       , col = coloursLegendENRICH
       , pch = 19
       , cex = 1.5)
dev.off()
 
jpeg("TaxaSummaries_ExNExN_Reduced.jpeg", width = 2000, height = 2000)
par(fig = c(0,01,0.5,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = colourLegendReduced
        , las = 2
        , border = NA
        , names.arg = sampleLab
        , space = c(0,0,0,0,0,1,0,0,0,1,0,0,0,0,2,0,0,0)
        , ylab = "Relative Abundance"
)
title(xlab = "Treatment", line = 6)
par(fig = c(0,1,0,0.5), mar = c(0,0,0,0), new = TRUE)
legend("center"
       , legend = OTUTAXONREDUCE
       , col = coloursLegendREDUCE
       , pch = 19
       , cex = 1.5)
dev.off()














# ORRR


#ExNWater
sampleLab <- c(rep("Meristem-Nereo",4), rep("Nereocystis",5), rep("Mastocarpus",5), rep("Nereo+Mast",5), rep("Water Only", 5))


# Get colour labels
OTUIDsENRICH <- rownames(colourLegendEnriched)[!is.na(colourLegendEnriched)]
OTUTAXONENRICH <- sapply(OTUIDs, function(x) {
  taxaID[grep(paste0("^",x,'$'), rownames(taxaID)),]
})
coloursLegendENRICH <- colourLegendEnriched[!is.na(colourLegendEnriched)]

OTUIDsREDUCE <- rownames(colourLegendReduced)[!is.na(colourLegendReduced)]
OTUTAXONREDUCE <- sapply(OTUIDs, function(x) {
  taxaID[grep(paste0("^",x,'$'), rownames(taxaID)),]
})
coloursLegendREDUCE <- colourLegendReduced[!is.na(colourLegendReduced)]



jpeg("TaxaSummaries_ExNWater_Enriched.jpeg", width = 2000, height = 2000)
par(fig = c(0,01,0.5,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = colourLegendEnriched
        , las = 2
        , border = NA
        , names.arg = sampleLab
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0)
        , ylab = "Relative Abundance"
)
title(xlab = "Treatment", line = 6)
par(fig = c(0,1,0,0.5), mar = c(0,0,0,0), new = TRUE)
legend("center"
       , legend = OTUTAXONENRICH
       , col = coloursLegendENRICH
       , pch = 19
       , cex = 1.5)
dev.off()

jpeg("TaxaSummaries_ExNWater_Reduced.jpeg", width = 2000, height = 2000)
par(fig = c(0,1,0.5,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = colourLegendReduced
        , las = 2
        , border = NA
        , names.arg = sampleLab
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0)
        , ylab = "Relative Abundance"
)
title(xlab = "Treatment", line = 6)
par(fig = c(0,1,0,0.5), mar = c(0,0,0,0), new = TRUE)
legend("center"
       , legend = OTUTAXONREDUCE
       , col = coloursLegendREDUCE
       , pch = 19
       , cex = 1.5)
dev.off()




###### PLOT ALL TAXASUMMARIES ######
jpeg("TaxaSummaries_ExNExN_all.jpeg", width = 1000, height = 1000)
par(fig = c(0,1,0,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = coloursRandom
        , names.arg = sampleLab
        , space = c(0,0,0,0,0,1,0,0,0,1,0,0,0,0,2,0,0,0)
        , las = 2
        , ylab = "Relative Abundance"
        , main = "Taxa Summaries-Meristem swabs "
        , border = FALSE
        )
title(xlab = "Treatment", line = 6)
dev.off()

jpeg("TaxaSummaries_ExNWater_all.jpeg", width = 1000, height = 1000)
par(fig = c(0,1,0,1), mar = c(8,4,4,1))
barplot(filteredAbundanceData
        , col = coloursRandom
        , names.arg = sampleLab
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0)
        , las = 2
        , ylab = "Relative Abundance"
        , main = "Taxa Summaries-Meristem swabs "
        , border = FALSE
)
title(xlab = "Treatment", line = 6)
dev.off()

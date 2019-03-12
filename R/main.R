# main.R
#
# Purpose: To automatically perform analysis on GEO experiments from JSON structured data
# Version: 1.0
# Date:    2019 03 12
# Author:  Yoonsik Park <yoonsik.park@mail.utoronto.ca>
#
# Input: JSON file
# Output: 
# Dependencies:
#
# Version history:
#
# References Used:
#         Starter Code from: https://github.com/hyginn/zu/blob/master/inst/scripts/expressionAnalysisSampleCode.R
# Notes:  R Code adapted from RPR-GEO2R  ABC learning unit.
#
#         Quantile normalization needs to cite:
#           Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
#           A Comparison of Normalization Methods for High Density
#           Oligonucleotide Array Data Based on Bias and Variance.
#           Bioinformatics 19(2) ,pp 185-193.
#
#
# ==============================================================================

# ====  PARAMETERS  ============================================================

# json filename
inputFile <- "sample.json"

# output filename for RData
outputFile <- "sampletable.RData"

# log file
logFile <- "sample.log"


# ====  PACKAGES  ==============================================================

if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}

if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}

# for json parsing
if (! require(rjson, quietly=TRUE)) {
  install.packages("rjson")
  library(rjson)
}


# for quantile normalization ...
if (! require(preprocessCore, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("preprocessCore")
  library(preprocessCore)
}

# JSON Description should be unique

# === Get HUGO symbols

load("./inst/extdata/HUGOsymbols.RData")
load("inst/extdata/synMap.RData")

convertToMatrix <- function(gse) {
  # We use the code from the GEO Vignette to recreate the ExpressionSet
  
  gsmlist = GSMList(gse)
  probesets <- Table(GPLList(gse)[[1]])$ID
  data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
  {tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
  }))
  data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
  data.matrix <- log2(data.matrix)
  data.matrix[1:5,]
  rownames(data.matrix) <- probesets
  colnames(data.matrix) <- names(gsmlist)
  pdata <- data.frame(samples=names(gsmlist))
  rownames(pdata) <- names(gsmlist)
  pheno <- as(pdata,"AnnotatedDataFrame")
  eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
  
  return(eset2)
}

getAveragedHGNCTable <- function(geoSeriesObj, geoSeriesObj2, description, controlSamples, experimentalSamples, isMicroArray=TRUE) {
  
  # Check if all Sample Names have been found in the GEO object
  
  allSampleNames <- names(GSMList(geoSeriesObj))
  for (sample in controlSamples) {
    if (!sample %in% allSampleNames) {
      warning(sample)
      stop("Sample from controlSamples not found!")
    }
  }
  for (sample in experimentalSamples) {
    if (!sample %in% allSampleNames) {
      warning(sample)
      stop("Sample from experimentalSamples not found!")
    }
  }

  eset <- convertToMatrix(geoSeriesObj)
  exprData <- exprs(eset)
  controls <- exprData[,unlist(controlSamples)]
  controls <- matrix(rowMeans(controls), ncol = 1)
  row.names(controls) <- row.names(exprData)
  colnames(controls) <- c("controls")
  
  experiments <- exprData[,unlist(experimentalSamples)]
  experiments <- matrix(rowMeans(experiments), ncol = 1)
  row.names(experiments) <- row.names(exprData)
  colnames(experiments) <- c("experiments")
  #combined <- cbind(controls, experiments)
  myAnnot <- NA
  if ("Gene symbol" %in% colnames(geoSeriesObj2@featureData@data)){ 
    myAnnot <- geoSeriesObj2@featureData@data[ , c("ID", "Gene symbol")]
  } else {
    myAnnot <- geoSeriesObj2@featureData@data[ , c("ID", "Symbol")]
  }
  
  colnames(myAnnot) <- c("probeIDs", "symbols") 
  myAnnot$probeIDs <- as.character(myAnnot$probeIDs)  # convert to character
  myAnnot$symbols  <- as.character(myAnnot$symbols)
  any(is.na(myAnnot$symbols))       # FALSE

  retGene <- matrix(nrow = length(HUGOsymbols))
  row.names(retGene) <- HUGOsymbols
  colnames(retGene) <- c("logRatioChange")
  for (gene in HUGOsymbols) {
    if (gene %in% myAnnot$symbol){ 
      sel <- which(myAnnot$symbol == gene)
      if (length(sel) > 1) {
        retGene[gene,1] <- log(mean(experiments[sel])/mean(controls[sel]))
      } else {
        retGene[gene,1] <- log(experiments[sel]/controls[sel])
      }
      
      }
  }
  return(retGene)
}

parseJson <- function(filename) {
  result <- fromJSON(file = filename)
  cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
  cat(": JSON Successfully Parsed\n", file = logFile, append=TRUE)
  first <- TRUE
  retData <- matrix(nrow = length(HUGOsymbols))
  row.names(retData) <- HUGOsymbols
  oldgName <- ""
  geoObj <- NA
  geoObj2 <- NA
  cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
  cat(paste(":",length(result)), file = logFile, append=TRUE)
  cat(" Experiments to be Analyzed\n", file = logFile, append=TRUE)
  counter <- 1
  for (experiment in result) {
    gName <- experiment$geoDataset
    d  <- experiment$description
    c  <-  as.list(experiment$controlSamples)
    e  <-  as.list(experiment$experimentalSamples)
    
    cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
    cat(paste(": Analyzing Experiment ", counter,"...", d), file = logFile, append=TRUE)
    cat("\n", file = logFile, append=TRUE)
    
    # Quick check to stop reloading the GEO data file unnecessarily
    if (oldgName != gName) {
      # We set GSEMatrix equal to FALSE because we need the GSM IDs
      geoObj <- getGEO(gName, GSEMatrix=FALSE, AnnotGPL=TRUE)
      # For some reason GSEMatrix=FALSE doesn't include the feature data
      geoObj2 <- getGEO(gName, GSEMatrix=TRUE, AnnotGPL=TRUE)[[1]]
    }
    oldgName <- gName

    newCol <- getAveragedHGNCTable(geoObj, geoObj2, d, c, e)
    colnames(newCol) <- d
    if (first) {
      first <- FALSE
      retData <- newCol
    } else {
      retData <- cbind(retData, newCol)
    }
    counter <- counter + 1
  }
  return(retData)
}

library("rjson")

# g <- "GSE35330"
# d  <- "Experiment 1"
# c  <-  list("GSM866121", "GSM866122", "GSM866124", "GSM866126")
# e  <-  list("GSM866111", "GSM866116", "GSM866118", "GSM866119")
# new <- getAveragedHGNCTable(g, d, c, e)

cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
cat(": Started JSON Parsing\n", file = logFile, append=TRUE)
ret <- parseJson(inputFile)

# Remove all rows with a NA
retNoNA <- ret[complete.cases(ret), ]

cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
cat(paste(": ", nrow(retNoNA), " out of ", nrow(ret), " HUGO genes were complete across all experiments\n"), file = logFile, append=TRUE)


afterQuantile <- normalize.quantiles(retNoNA, copy = FALSE)

# Finally save the RData file

save(afterQuantile, file = outputFile)
cat(format(Sys.time(), "%a %b %d %X %Y"), file = logFile, append=TRUE)
cat(paste(": RData Succesfully Writtten to ", outputFile,"\n"), file = logFile, append=TRUE)

# [END]

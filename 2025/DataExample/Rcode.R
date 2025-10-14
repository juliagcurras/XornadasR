###############################################################################-
###################                                       #####################-
###################           RAW DATA ANALYSIS           #####################-
###################                                       #####################-
###############################################################################-


# Julia G Currás  - 2025/10/13

# Setup ####
setwd("H:/Mi unidad/Proyectos/Bioinformatica/PipelineOmics/DataExample")
library(dplyr)
library(ggplot2)
set.seed(9396)

doTestT <- function(df, dfGrupos, g1, g2) {
  
  # g1: control
  df <- as.data.frame(df)
  dfGrupos <- as.data.frame(dfGrupos)
  
  # Selecting columns from dataframe of each group
  samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, "Samples"])
  samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, "Samples"])
  df <- df[, c(samplesG1, samplesG2)]
  
  # log2FC: group2 (case) - group1 (control)
  df$FC  <- apply(df, 1, function(x) {
    mean(x[samplesG2], na.rm =T)/mean(x[samplesG1], na.rm = T)
  }
  )
  df$logFC <- log(df$FC, base = 2)
  
  # Fifference of means estimation
  df$DiffExpr  <- apply(df, 1, function(x) {
    (mean(x[samplesG2], na.rm =T) - mean(x[samplesG1], na.rm = T))
  }
  )
  
  # T-test with equal variance
  df[, c( "t", "P.Val")] <- t(apply(df, 1, function(x) {
    # When there isn't enough observations, returning NA
    evalG1NA <- sum(is.na(x[samplesG1])) #/length(x[samplesG1]) 
    evalG2NA <- sum(is.na(x[samplesG2]))
    if (((length(x[samplesG1])-evalG1NA) < 3) | 
        ((length(x[samplesG2])-evalG2NA) < 3 ) #|
    ) {
      return(c(NA, NA))
    }
    # Otherwise, test
    testDone <- try(stats::t.test(x[samplesG2], x[samplesG1],
                                  alternative = "two.sided",
                                  var.equal = TRUE), T)
    if (!(class(testDone) == "try-error")){ # if some error happend, then NA
      res <- testDone
      return(c(res$statistic, res$p.value))
    } else {
      return(c(NA, NA)) 
    }
  }
  ))
  
  #Benjamini-Hochberg correction for multiple testing
  df$adj.P.Val <- stats::p.adjust(df$P.Val, method = "BH")
  df <- df[, c("logFC", "FC", "P.Val", "adj.P.Val", "DiffExpr")]
  
  # Returning final vector
  return(df)
}


removeDiscordancesReplicates <- function(dfQ, nReps = 2){
  posFinal <- nReps
  dfResult <- NULL
  for (i in seq(1, ncol(dfQ)-1, nReps)){ # for each sample
    if (i != 1){
      posFinal <- posFinal + nReps
    } 
    dfAux <- dfQ[, i:posFinal] # its respective columns are selected
    
    resultProt <- apply(dfAux, 1, function(fila){ # and then only one value is established
      totalNA <- sum(is.na(fila)) 
      if (totalNA >=3){ # NA if there are 3 or 4 missing values
        return(NA)
      } else {
        return(mean(fila, na.rm = T)) # mean(replicates) if there are 0, 1 or 2 missing values
      }
    })
    
    # Finally, info from different samples is bound
    if (is.null(dfResult)){
      dfResult <- data.frame(resultProt)
    } else {
      dfResult <- cbind(dfResult, resultProt)
    }
  }
  colnames(dfResult) <- colnames(dfQ)[(seq(nReps, ncol(dfQ), nReps))]
  
  # Final object!
  return(dfResult)
}


## Reanalysis ####
# Loading data from TS5 #
tab5 <- xlsx::read.xlsx(file = "377_100_PubMed_table5.XLSX", sheetIndex = 1)
colnames(tab5) <- tab5[1,] 
tab5 <- tab5[-1,]

# Selecting study groups #
data <- tab5 %>% 
  dplyr::select(Accession_id, Description, starts_with("Abundance_DM"), starts_with("Abundance_NC"))

# Keeping only abundance information #
dfQuant <- data[, -c(1,2)]
rownames(dfQuant) <- data$Accession_id
colnames(dfQuant) <- gsub(colnames(dfQuant), pattern = "Abundance_", replacement = "")
dfQuant[dfQuant == 0] <- NA

# Design matrix #
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = gsub(colnames(dfQuant), pattern = "[0-9]", replacement = "")
)


out <- list(
  dataMatrix = dfQuant,
  designMatrix = dfGrupos
)
saveRDS(out, file = "rawDataMatrix.rds")




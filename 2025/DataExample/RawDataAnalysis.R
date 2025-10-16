###############################################################################-
###################                                       #####################-
###################           RAW DATA ANALYSIS           #####################-
###################                                       #####################-
###############################################################################-


# Julia G Currás  - 2025/10/13

# Setup ####
setwd("C:/Users/julia/Documents/GitHub/XornadasR/2025/DataExample")
library(dplyr)
library(tidyr)
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

cyclicLoessNorm <- function (log2Matrix) {
  log2Matrix <- as.matrix(log2Matrix)
  normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method = "fast")
  colnames(normMatrix) <- colnames(log2Matrix)
  # Remove is.infinite data
  normMatrix[sapply(normMatrix, is.infinite)] <- NA
  normMatrix
}





## Reanalysis ####
# Loading data from TS5 #
data <- read.table(file = "proteinGroups.txt", 
                   header = T, sep = "\t")
df <- data %>% dplyr::select(Protein.IDs, Protein.names, Gene.names, 
                             starts_with("LFQ.intensity."), Reverse)
# Set 0 to NA 
df[df == 0] <- NA

# filter
df <- df %>%
  dplyr::filter(Reverse != "+") %>%
  dplyr::select(-Reverse) 
colnames(df) <- gsub(x = colnames(df), pattern = "LFQ.intensity.", replacement = "")

# selecting samples from group 2 (T2DM) and 4 (control) #
dfQuant <- df %>% dplyr::select(Group4_rep1:Group4_rep8, Group2_rep1:Group2_rep8) # 2vs4
colnames(dfQuant) <- gsub(pattern = "Group2", replacement = "T2DM", x = colnames(dfQuant))
colnames(dfQuant) <- gsub(pattern = "Group4", replacement = "Control", x = colnames(dfQuant))

# Building design matrix #
rownames(dfQuant) <- df$Protein.IDs
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = substr(colnames(dfQuant), start = 1, stop = nchar(colnames(dfQuant)) - 5)
)
table(dfGrupos$Groups)

# Excluding proteins with more than 2 missing values in any group #
g1 <- "Control"
g2 <- "T2DM"
naLim <- 4/8
samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, "Samples"])
samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, "Samples"])
indOkProts <- apply(dfQuant, 1, function(x) {
  evalG1NA <- sum(is.na(x[samplesG1]))/length(x[samplesG1]) 
  evalG2NA <- sum(is.na(x[samplesG2]))/length(x[samplesG2])
  return(!any(evalG1NA >= naLim, evalG2NA >= naLim)) # returning only the good ones
})
dfQuant <- as.data.frame(dfQuant)
dfQuant1 <- dfQuant[which(indOkProts),]

# Log transformation #
df <- as.data.frame(log(dfQuant1, base = 2))

cores <- grDevices::colorRampPalette(colors = c("#7CA9B2", "#23373B","#0A0F10"))(16)
# Representacións gráficas
  # Sen normalizar #
dfQuantLog <- log(dfQuant, base = 2)
datos_largo <- dfQuantLog %>%
  pivot_longer(
    cols = colnames(df), # Selecciona todas las columnas que empiezan con 'Dato_'
    names_to = "Mostras", # Nombre para la nueva columna de los identificadores
    values_to = "Cantidade" # Nombre para la nueva columna con los valores
  )

pLog <- ggplot(datos_largo, aes(x = Mostras, y = Cantidade, fill = Mostras)) +
  geom_boxplot(
    outlier.shape = 1, # Forma de los puntos atípicos (pequeños círculos)
    outlier.size = 0.5, # Tamaño de los puntos atípicos
    alpha =0.7 # Color del borde de las cajas
  ) +
  ggplot2::scale_fill_manual(values = cores) +
  ggplot2::scale_color_manual(values = cores) +
  labs(
    title = "Sen normalizar",
    # x = "Mostras",
    y = "Cantidade de proteínas"
  ) +
  ylim(c(20,40))+
  theme_minimal(base_size = 12) +
  theme(
    axis.line = ggplot2::element_line(linewidth = 0.5,colour = "black"), 
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5), legend.position = "none"
  )

# pLog <- Biostatech::plotBoxMultivar(base = dfQuantLog, varResumen = colnames(dfQuantLog), interact = T,,
#                                     tituloX = "Sen normalizar", tituloY ="Cantidade de proteínas", #tituloX = "Mostras",
#                             color = grDevices::colorRampPalette(colors = c("#7CA9B2", "#23373B","#0A0F10"))(43))$grafico

# Normalizando #
df <- cyclicLoessNorm(df)
df <- as.data.frame(df)
datos_largo <- df %>%
  pivot_longer(
    cols = colnames(df), # Selecciona todas las columnas que empiezan con 'Dato_'
    names_to = "Mostras", # Nombre para la nueva columna de los identificadores
    values_to = "Cantidade" # Nombre para la nueva columna con los valores
  )

pNorm <- ggplot(datos_largo, aes(x = Mostras, y = Cantidade, fill = Mostras)) +
  geom_boxplot(
    outlier.shape = 1, # Forma de los puntos atípicos (pequeños círculos)
    outlier.size = 0.5, # Tamaño de los puntos atípicos
    alpha =0.7  # Color del borde de las cajas
  ) +
  ggplot2::scale_fill_manual(values = cores) +
  ggplot2::scale_color_manual(values = cores) +
  labs(
    title = "Normalizados - CyclicLoess",
    # x = "Mostras",
    y = "Cantidade de proteínas"
  ) +
  ylim(c(20,40))+
  theme_minimal(base_size = 12) +
  theme(
    axis.line = ggplot2::element_line(linewidth = 0.5,colour = "black"), 
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5), legend.position = "none"
  )
# pNorm <- Biostatech::plotBoxMultivar(
#   base = df,
#   varResumen =colnames(df), interact = T,
#   tituloX = "Normalizado - CyclicLoess", tituloY ="Cantidade de proteínas", #tituloX = "Mostras",
#   color = grDevices::colorRampPalette(colors = c("#7CA9B2", "#23373B", "#0A0F10"))(43))$grafico
# outGraphs <- list(
#   Raw = pLog,
#   Norm = pNorm
# )
# saveRDS(outGraphs, file = "preprocessGraphs.rds")

combined <- ggpubr::ggarrange(pLog, pNorm, ncol = 1)
# combined <- subplot(pLog, pNorm, nrows = 2, shareY = F, shareX = F, titleY = T, titleX = T, margin = 0.05)
saveRDS(combined, file = "preprocessGraphs_.rds")

# Saving
out <- list(
  dataMatrix = dfQuant,
  normDataMatrix = df,
  designMatrix = dfGrupos
)
saveRDS(out, file = "rawDataMatrix.rds")



#.--------------------------------------------------------------------------####s
out <- readRDS(file = "rawDataMatrix.rds")
df <- out$normDataMatrix
infoName <- data %>% dplyr::select(Protein.IDs, Gene.names)
# Expresión diferencial ####
tabRes <- doTestT(df = df, dfGrupos = dfGrupos, g1 = g1, g2 = g2)
tabRes$ID <- rownames(tabRes) 
tabRes$Description <- data[which(data$Accession_id %in% tabRes$ID), "Description"]

dfRes <- merge(tabRes, infoName, by.x = "ID", by.y = "Protein.IDs")
dfRes <- dfRes %>% dplyr::select(ID, Gene.names, everything(), -DiffExpr)
colnames(dfRes) <- c("ID", "Proteina", "logFC", "FC", "P-valor", "P-valor ajustado")
dfRes$logFC <- dfRes$logFC*20
dfRes$FC <- dfRes$FC*20
rownames(dfRes) <- dfRes$ID
dfRes <- dfRes %>% arrange(`P-valor`)

saveRDS(dfRes, file = "resultsDE.rds")



#.--------------------------------------------------------------------------####
out <- readRDS(file = "rawDataMatrix.rds")
dfRes <- readRDS(file = "resultsDE.rds")
df <- out$normDataMatrix
a <- 5
signProts <- dfRes %>% dplyr::filter(`P-valor`<0.05) %>% pull(ID)
df <- df[signProts,]
df[c(1, 4:7, 11, 16, 17), 1:8] <- df[c(1, 4:7, 11, 16, 17), 1:8] + a
df[c(1, 4:7, 11, 16, 17), 9:16] <- df[c(1, 4:7, 11, 16, 17), 9:16] - a
df[c(2:3, 8:10, 12:15), 1:8] <- df[c(2:3, 8:10, 12:15), 1:8] -a
df[c(2:3, 8:10, 12:15), 9:16] <- df[c(2:3, 8:10, 12:15), 9:16] +a

df <- as.data.frame(t(df))
df$Grupos <- dfGrupos$Groups
colnames(df) <- sapply(strsplit(colnames(df), split = ";"), "[[", 1)

saveRDS(df, file = "dfClust.rds")

# CLUSTERING ####
heatmaply::heatmaply(df, scale = "row", k_row = 2, 
                     color = c("green", "black", "red"))

pcaRes <- pcaMethods::pca(object = as.matrix(df[,-ncol(df)]), 
                       method = "nipals", nPcs = 2, 
                       scale = "none", center = F)

# Adding group info to plot
dfPlotPCA <- merge(df, pcaMethods::scores(pcaRes), by=0)

#### Clustering por grupos 

ggplot(dfPlotPCA, aes(PC1, PC2, colour=Grupos)) +
  geom_point() +
  stat_ellipse() + 
  theme_minimal() +
  labs(colour = "Grupos") +
  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
                 axis.ticks = ggplot2::element_line(linewidth = 0.5, 
                                                    colour = "black"))



dfPCA <- t(na.omit(t(df[,-ncol(df)])))
res.pca <- prcomp(dfPCA, scale = FALSE)
# PCA main figure
a <- 5
factoextra::fviz_pca_ind(res.pca,
             habillage = dfGrupos$Groups,
             palette = c("#9EBFC6", "#23373B"),
             addEllipses = TRUE,
             label = "none",
             title = "PCA") + 
  theme_minimal() +
  theme(text = element_text(size = 16-a),
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        title = element_text(size = 20-a),
        axis.title = element_text(size = 18-a),
        axis.text = element_text(size = 16-a),
        legend.text = element_text(size = 16-a),
        legend.title = element_text(size = 18-a))



#.--------------------------------------------------------------------------####
out <- readRDS(file = "rawDataMatrix.rds")
dfRes <- readRDS(file = "resultsDE.rds")
# Análisis funcional ####
library(org.Hs.eg.db)

protsID <- dfRes %>% filter(`P-valor`<0.2) %>% pull(ID)
resultGO <- clusterProfiler::enrichGO(
  gene = protsID,
  OrgDb = "org.Hs.eg.db",
  keyType = "UNIPROT",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  # universe = protsIDUniverse,
  qvalueCutoff = 0.2,
  minGSSize = 3,
  maxGSSize = 800,
  readable = TRUE,
  pool = FALSE
)
resultGO@result %>% View
saveRDS(resultGO, file = "GO_ORA.rds")
resultGO <- readRDS(file = "GO_ORA.rds")

library(enrichplot)

resultGO@result <- resultGO@result %>% filter(ONTOLOGY == "BP")
enrichplot::dotplot(resultGO, showCategory = 10)
enrichplot::cnetplot(resultGO, color_gene = "darkolivegreen4", 
                     color_category = "purple", circular=F, 
                     showCategory = 10)






# Alternativa heatmap ####
library(ggplot2)
df <- readRDS(file = "dfClust.rds")
df <- df[,-ncol(df)]
df <- na.omit(df)
# Matrix of common protein to long df
df_long <- reshape2::melt(as.matrix(df))
colnames(df_long) <- c("Sample", "Protein", "Abundance")

# Using the scaled importance to find the optimal order of the proteins based 
# on the similarity between methods (seriation package)
score_matrix_scaled <- reshape2::acast(df_long, Protein ~ Sample, 
                                       value.var = "Abundance")
abc <- rownames(score_matrix_scaled) 
score_matrix_scaled <- apply(score_matrix_scaled, 2, as.numeric)
rownames(score_matrix_scaled)  <- abc
order_rows <- seriation::get_order(seriation::seriate(score_matrix_scaled, method = "PCA"))
# Reordering 
ordered_proteins <- rownames(score_matrix_scaled)[order_rows]
df_long$Protein <- factor(df_long$Protein, levels = ordered_proteins)

# Creating heatmap with reorder rows. Values on cells are non-scaled importance. 
p1 <- ggplot2::ggplot(df_long, aes(x=Protein, y=Sample, fill=Abundance)) +
  ggplot2::geom_tile(color="white") +
  # ggplot2::geom_text(aes(label=round(Abundance, 2)), size=15) +
  ggplot2::scale_fill_gradient(
    low = "#A6CAEC", 
    high = "#08306B") +
  ggplot2::theme_minimal() +
  # ggplot2::labs(
  #   title="",
  #   fill="Abundance") +
  xlab("") + ylab("") +
  ggplot2::theme(
    axis.text.y = element_text(face = "italic", size = 40),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 45), 
    legend.position = "bottom",
    legend.key.size = unit(4.3, "cm"),
    legend.title = element_text(size = 40),    
    legend.text = element_text(size = 40) 
  )
p1

















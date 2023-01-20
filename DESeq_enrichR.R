library('tidyverse')
library('tximport')
library('ensembldb')
library('EnsDb.Hsapiens.v86')
library('RUVSeq')
library('data.table')
library('edgeR')
library('matrixStats')
library('cowplot')
library('EDASeq')
library('ggforce')
library('gt')
library('plotly')
library('NOISeq')
library('biomaRt')
library('GenomicFeatures')
library('useful')
library('baySeq')
library('limma')
library('DESeq2')
library('TCC')
library('snow')
library('DEGseq')
library('enrichR')
library('biomaRt')
library('foreign')
library('xtable')
library('stargazer')
library('org.Hs.eg.db')

setwd("~/Escritorio/RESULTS/CART14")
metadata <- as.data.frame(read_csv('metadata_deseq.csv'))
#metadata <- metadata[-c(3,7),]
path <- list.files(".", pattern = "*_ReadsPerGene.out.tab")
files_names <- paste0(sapply(strsplit(path, '_'), `[[`, 1))
group <- metadata$prog
# Load read per genes files
countData_LT <- data.frame(fread(path[1]))[c(1,2)]
for(i in 2:length(path)) {
  countData_LT <- cbind(countData_LT, data.frame(fread(path[i]))[c(2)])
}
countData_LT <- countData_LT[c(5:nrow(countData_LT)),]
countData_LT <- data.frame(countData_LT, row.names = c(1))
setnames(countData_LT, files_names)

dge_LT.raw <- DGEList(counts = countData_LT, group = group)


log2.cpm <- cpm(dge_LT.raw, log=TRUE)
log2.cpm.df <- as.tibble(log2.cpm, rownames = "ensID")
log2.cpm.df
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = files_names[1]:files_names[i],
                                  names_to = 'samples',
                                  values_to = 'expression')
nonnon <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = 'median',
               geom = 'point',
               shape = 95,
               size = 10,
               color = 'black',
               show.legend = FALSE) +
  labs(y='log2 expression', x = 'sample',
       title = 'log2 Counts per Million (CPM)',
       subtitle = 'unfiltered, non-normalized') +
  theme(axis.text.x = element_text(angle=40, vjust = 0.5, hjust = 0.5)) ##mover antes del filtrado


## Filtrado
keep <- rowSums(dge_LT.raw$counts>0)>=i
dge.filt <- dge_LT.raw[keep,]
dim(dge.filt)
dim(dge_LT.raw)

log2.cpm.filt <- cpm(dge.filt, log=TRUE)
log2.cpm.filt.df <- as.tibble(log2.cpm.filt, rownames = "ensID")
log2.cpm.filt.df
log2.cpm.filt.df.pivot <- pivot_longer(log2.cpm.filt.df,
                                       cols = files_names[1]:files_names[i],
                                       names_to = 'samples',
                                       values_to = 'expression')
yesnon <- ggplot(log2.cpm.filt.df.pivot) +
  aes(x=samples, y=expression, fill=samples) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = 'median',
               geom = 'point',
               shape = 95,
               size = 10,
               color = 'black',
               show.legend = FALSE) +
  labs(y='log2 expression', x = 'sample',
       title = 'log2 Counts per Million (CPM)',
       subtitle = 'filtered, non-normalized') +
  theme(axis.text.x = element_text(angle=40, vjust = 0.5, hjust = 0.5))


## Normalización TMM
dge.filt.norm <- calcNormFactors(dge.filt, method = 'TMM')
log2.cpm.filt.norm <- cpm(dge.filt.norm, log=TRUE)
log2.cpm.filt.norm.df <- as.tibble(log2.cpm.filt.norm, rownames = "ensID")
log2.cpm.filt.norm.df
log2.cpm.filt.norm.df.pivot <- pivot_longer(log2.cpm.filt.norm.df,
                                            cols = files_names[1]:files_names[i],
                                            names_to = 'samples',
                                            values_to = 'expression')
yesyes <- ggplot(log2.cpm.filt.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = 'median',
               geom = 'point',
               shape = 95,
               size = 10,
               color = 'black',
               show.legend = FALSE) +
  labs(y='log2 expression', x = 'sample',
       title = 'log2 Counts per Million (CPM)',
       subtitle = 'filtered, normalized') +
  theme(axis.text.x = element_text(angle=40, vjust = 0.5, hjust = 0.5))
plot_grid(nonnon, yesnon, yesyes, nrow = 2, ncol = 2,	label = 'AUTO')


###  DESeq2. SI FC POSITIVO <-> OVEREXPRESSED EN PROGRESION O RECAIDA
dds <- DESeqDataSetFromMatrix(countData = dge.filt, colData = metadata, design = ~ prog)
keep <- rowSums(counts(dds)>5)>=5 ##Reducir por si infraexpresion en recaidas.
dds.filt <- dds[keep,]
dds.filt <- DESeq(dds.filt)
res <- results(dds.filt)
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.1)
summary(res, res$padj<0.1)
write.csv(as.data.frame(resSig), './prog_res_deseq.csv')

#PCA
vsd <- vst(dds.filt, blind=FALSE)
mat <- assay(vsd)
mm <- model.matrix(~prog, colData(vsd))
mat <- limma::removeBatchEffect(mat, design=mm)
assay(vsd) <- mat

pcaData <- plotPCA(vsd, intgroup=c("prog"), returnData=TRUE, ntop = 1400)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=prog)) +
  geom_point(size=3) +
  #geom_label(aes(label = metadata$sample_id)) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(title='PCA plot', subtitle = "CAR-T cells analysis", colour = "Progression") +
  coord_equal(1.8) +
  scale_color_manual(name = 'Progression', values = c('lightskyblue3', 'salmon'), limits = c('no', 'yes'))

png(file='./pca.png')
plot(pca)
dev.off()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ENSGID <- row.names(resSig)
mygenes <- getBM(attributes = 'external_gene_name', filters = 'ensembl_gene_id', values = ENSGID,
                 mart = ensembl)
write.csv(mygenes, './genes.csv')
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
dbs <- c("Transcription_Factor_PPIs", "BioCarta_2016",
         "Human_Gene_Atlas", "Pfam_InterPro_Domains", "Reactome_2016",
         "NCI-Nature_2016", "Panther_2016", "GO_Molecular_Function_2018", 
         "GO_Biological_Process_2018", "TRANSFAC_and_JASPAR_PWMs")
if (websiteLive) {
  enriched <- enrichr(mygenes[,1], dbs)
}



if (websiteLive) enriched[["Transcription_Factor_PPIs"]] ### Los receptores de estrogenos ESR1-2 NO son fundamentales para el desarrollo T, sin embargo participa en la inflamaciñon mediada por célula T (en concreto Th1 y Th17, su underexpression reduce IFNgamma y reduce la act de células T). Asimismo ESR1 underxpression mejora la supervivencia celular al reducir la apoptosis mediadada por la activación de ceulas T. SU UNDEREXPRESSION REDUCE LA PROLIFERACION Y EXPANSIÓN de T cel
### ILF2-3 (Interleukin enhancer-binding factor 2) participan en procesos celulares como la replicación del ADN, reparación, estabilización del mRNA, inhibicion de transcripción y sintesis de miRNA. Están altamente relacionadas con las vias MAPK y PI3K. Sin ILF2 disminuye proliferación y aumenta apoptosis (overexpression al revés)
### TCF3 (EA2) mantiene el estado de doble positiva a las células T (estado de desarrollo) mientras que la presencia de RUNX hace que se desarrollen a CD8. Una overexpression de RUNX3 se ve en rec (por tanto +CD8 en prog que no prog)
### SMAD2 regula la síntesis de TGFB, regulando la proliferation, apoptosis y diferenciacion celular. En concreto aumenta las Treg que suponen que en overexpresion las cñelulas CART supriman su respuesta.
### NFKB aumenta supervivencia y la producción de citoquinas proinflamatorias.
### Downregulación de SMAD4 empeora la función T en autoimunidad y antitumoral.
if (websiteLive) enriched[["BioCarta_2016"]] ### Las paths que han sido significativas se relacionan con la activación de células T y su señalización por TCR. Asimismo, otros genes hacen overlaping con la apoptosis mediada MEF2D.
if (websiteLive) enriched[["Human_Gene_Atlas"]] ### Son células CD4 y CD8. 
if (websiteLive) enriched[["Pfam_InterPro_Domains"]] ### De aquí se observa que muchas proteínas de las overexpressed en recaida modulan union a RNA. 
if (websiteLive) enriched[["Reactome_2016"]] ### Se ven terminos significativos sobre la traducción, presentación de antigenos MHCI, señalización TCR, procesamiento y presentación antiǵenica, paths de señalización de TCR (incluyendo fosforilación TCRtheta), 
if (websiteLive) enriched[["NCI-Nature_2016"]] ### Significativa la señalización TCR en CD8 y CD4, así como señalización RhoA (act cél T) y TGFB
if (websiteLive) enriched[["Panther_2016"]] ### Significativo la activación de células T
if (websiteLive) enriched[["GO_Molecular_Function_2018"]] ### La mayor parte de los terminos hace referencia a la unión de prot-mRNA y a la ubiquitinación
if (websiteLive) enriched[["GO_Biological_Process_2018"]] ### Terminos de presentación de anígenos, ribosómicos, presentación antigénica, splicing, señalización de receptor de células T y reg. del sistema inmune.

##### HACER TABLAS
a <- as.data.table(enriched[2])[,c(1,2,4,9)]
gtab <- gt(as.data.table(enriched[9])[,c(1,2,4,9)], auto_align = TRUE)
gtab <- cols_label(gtab, Transcription_Factor_PPIs.Term = "Transcription Factor", 
                   Transcription_Factor_PPIs.Overlap = "Number of overlapping genes",
                   Transcription_Factor_PPIs.Adjusted.P.value = "Adjusted P-value", 
                   Transcription_Factor_PPIs.Genes = "Overlapping genes") %>% 
  tab_header(title = "Transcriptional Factors Protein-Protein Interactions") %>%
  tab_options(table.align = "center", heading.align = "center") %>%
  fmt_number(columns = 3,
             decimals = 2,
             suffixing = TRUE) %>%
  cols_align(align = "center", columns = everything())



library(clusterProfiler)
library(wordcloud)
library(org.Hs.eg.db)

original_gene_list <- resSig$log2FoldChange

# name the vector
names(original_gene_list) <- row.names(resSig)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(resSig, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- row.names(sig_genes_df)

# omit NA values
genes <- na.omit(resSig)


resSig$symbol <- mapIds(org.Hs.eg.db, keys=row.names(resSig), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
results_gene_name <- as.data.frame(resSig)
write.csv(as.data.frame(results_gene_name), './results_gene_name.csv')

result_input_GO <- results_gene_name[c(2,6,7)]
result_input_GO <- data.frame(result_input_GO$symbol, result_input_GO$padj)
result_input_GO <- na.omit(result_input_GO)
RA_processed_1 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "KEGG")
RA_processed_2 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "Reactome")
RA_processed_3 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "BioCarta")
RA_processed_4 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "GO-AII")
RA_processed_5 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "GO-BP")
RA_processed_6 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "GO-CC")
RA_processed_7 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "GO-MF")
RA_processed_8 <- run_pathfindR(result_input_GO, p_val_threshold = 0.1, gene_sets = "cell_markers")


enriched_PPI <- enriched[[1]][c(1,2,4,9)] 
overlap <- enriched[[1]][2]
enriched_PPI$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_PPI$Down_regulated <- NA
enriched_PPI[is.na(enriched_PPI)] <- ""
colnames(enriched_PPI) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_PPI$Up_regulated <- gsub('\\;', ', ', enriched_PPI$Up_regulated)
a <- enrichment_chart(enriched_PPI, top_terms = 26)
png(file='./PPI.png')
plot(a)
dev.off()

enriched_BioCarta <- enriched[[2]][c(1,2,4,9)] 
enriched_BioCarta$Term <- str_split_i(enriched_BioCarta$Term, pattern = " Homo", 1)
overlap <- enriched[[2]][2]
enriched_BioCarta$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_BioCarta$Down_regulated <- NA
enriched_BioCarta[is.na(enriched_BioCarta)] <- ""
colnames(enriched_BioCarta) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_BioCarta$Up_regulated <- gsub('\\;', ', ', enriched_BioCarta$Up_regulated)
b <- enrichment_chart(enriched_BioCarta, top_terms = 15) #
png(file='./BioCarta.png')
plot(b)
dev.off()

enriched_GA <- enriched[[3]][c(1,2,4,9)] 
enriched_GA$Down_regulated <- NA
overlap <- enriched[[3]][2]
enriched_GA$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_GA[is.na(enriched_GA)] <- ""
colnames(enriched_GA) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_GA$Up_regulated <- gsub('\\;', ', ', enriched_GA$Up_regulated)
c <- enrichment_chart(enriched_GA, top_terms = 4)
png(file='./GeneAtlas.png')
plot(c)
dev.off()

enriched_reac <- enriched[[5]][c(1,2,4,9)] 
enriched_reac$Term <- str_split_i(enriched_reac$Term, pattern = " Homo", 1)
overlap <- enriched[[5]][2]
enriched_reac$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_reac$Down_regulated <- NA
enriched_reac[is.na(enriched_reac)] <- ""
colnames(enriched_reac) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_reac$Up_regulated <- gsub('\\;', ', ', enriched_reac$Up_regulated)
Reactome <- enrichment_chart(enriched_reac, top_terms = 50) 
png(file='./reactome.png')
plot(Reactome)
dev.off()

enriched_nci <- enriched[[6]][c(1,2,4,9)] 
enriched_nci$Term <- str_split_i(enriched_nci$Term, pattern = " Homo", 1)
overlap <- enriched[[6]][2]
enriched_nci$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_nci$Down_regulated <- NA
enriched_nci[is.na(enriched_nci)] <- ""
colnames(enriched_nci) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_nci$Up_regulated <- gsub('\\;', ', ', enriched_nci$Up_regulated)
NCI <- enrichment_chart(enriched_nci, top_terms = 24) #
png(file='./NCI.png')
plot(NCI)
dev.off()


enriched_pan <- enriched[[7]][c(1,2,4,9)] 
enriched_pan$Term <- str_split_i(enriched_pan$Term, pattern = " Homo", 1)
enriched_pan$Down_regulated <- NA
overlap <- enriched[[7]][2]
enriched_pan$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_pan[is.na(enriched_pan)] <- ""
colnames(enriched_pan) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_pan$Up_regulated <- gsub('\\;', ', ', enriched_pan$Up_regulated)
panther <- enrichment_chart(enriched_pan, top_terms = 11) #
png(file='./panther.png')
plot(panther)
dev.off()

enriched_MF <- enriched[[8]][c(1,2,4,9)] 
enriched_MF$Down_regulated <- NA
overlap <- enriched[[8]][2]
enriched_MF$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_MF[is.na(enriched_MF)] <- ""
colnames(enriched_MF) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_MF$Up_regulated <- gsub('\\;', ', ', enriched_MF$Up_regulated)
GO_MF <- enrichment_chart(enriched_MF, top_terms = 18) #
png(file='./GO-MF.png')
plot(GO_MF)
dev.off()

enriched_BF <- enriched[[9]][c(1,2,4,9)] 
enriched_BF$Down_regulated <- NA
overlap <- enriched[[9]][2]
enriched_BF$Overlap <- sapply(overlap[,], function(x) eval(parse(text=x)))
enriched_BF[is.na(enriched_BF)] <- ""
colnames(enriched_BF) <- c("Term_Description", "Fold_Enrichment", "lowest_p", "Up_regulated", "Down_regulated")
enriched_BF$Up_regulated <- gsub('\\;', ', ', enriched_BF$Up_regulated)
GO_BF <- enrichment_chart(enriched_BF, top_terms = 50)
png(file='./GO-BF.png')
plot(GO_BF)
dev.off()

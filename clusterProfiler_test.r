# Introduction: Gene Set Enrichment Analysis (GSEA) using the clusterProfiler package
# For detailed information, refer to the official documentation: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# Check and install required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "pathview", "enrichplot"))
require(DOSE)

# Load necessary libraries
required_packages <- c("clusterProfiler", "enrichplot", "ggplot2", "ggnewscale", "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes", "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra")

# Load each package separately
for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
    }
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
    }
    library(package, character.only = TRUE)
}

# Install and load packages
install_load_packages <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
    }
    library(package, character.only = TRUE)
}
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(ggridges)
library(europepmc)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(ggthemes)
library(ggExtra)
library(ggforce)
library(ggalluvial)
library(lattice)
library(latticeExtra)

# Check if packages are installed, if not install them
lapply(required_packages, install_load_packages)

# Annotation setup: Select the organism and install the corresponding annotation package
organism <- "org.Mm.eg.db"  # Set the desired organism code
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# set working directory
setwd("S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA")
# Define Results dir S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA/Results
results_dir <- "S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA/Results"

# define cell types to compare
cell_types <- c("neuropil", "microglia")

# Prepare Input: Read and process data
# Automatically read in cells from cell types and adjust the .csv
df <- read.csv(paste0("S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA/Datasets/", paste(cell_types, collapse = "_"), ".csv"), header = TRUE)

colnames(df)[1] <- "gene_symbol"
original_gene_list <- df$log2fc
names(original_gene_list) <- df$gene_symbol
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none")

# Visualizations
# Dotplot
require(DOSE)
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
# add title to dotplot
p1 <- dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
p1 <- p1 + labs(title = paste("GSEA of", paste(cell_types, collapse = " over ")))
p1
# save dotplot to results dir
ggsave(paste0(results_dir, "/GSEAdotplot_", paste(cell_types, collapse = "_"), ".png"), p1, units = "cm", dpi = 300)

# Enrichment Plot Map
similarity_matrix <- pairwise_termsim(gse)
emapplot(similarity_matrix, showCategory = 10)

# Category Netplot
cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)

# Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

# GSEA Plot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# PubMed Trend of Enriched Terms
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion = FALSE)

# KEGG Gene Set Enrichment Analysis
# Convert gene IDs for gseKEGG function
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "UNIPROT", OrgDb = organism)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
df2 <- df[df$gene_symbol %in% dedup_ids$SYMBOL,]
df2$gene_symbol <- dedup_ids$UNIPROT
kegg_gene_list <- df2$log2fc
names(kegg_gene_list) <- df2$gene_symbol
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

# KEGG GSEA
kegg_organism <- "mmu"
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = kegg_organism,
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "uniprot")

# Visualizations for KEGG GSEA KEGG GSEA Enriched Pathways of
p2 <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split = ".sign") + facet_grid(. ~ .sign)
p2 <- p2 + labs(title = paste("KEGG GSEA Enriched Pathways of", paste(cell_types, collapse = " over ")))
p2
# remove dot contour lines
# save dotplot to results dir
ggsave(paste0(results_dir, "/KEGGpathway_", paste(cell_types, collapse = "_"), ".png"), p2, units = "cm", dpi = 300)

emapplot(kk2)
cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

# Pathview: Create enriched KEGG pathway plots
install.packages("pathview")
library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data = kegg_gene_list, pathway.id = "dme04130", species = kegg_organism)
# Produce a different plot (PDF)
dme <- pathview(gene.data = kegg_gene_list, pathway.id = "dme04130", species = kegg_organism, kegg.native = FALSE)

# Display the native KEGG plot
knitr::include_graphics("dme04130.pathview.png")

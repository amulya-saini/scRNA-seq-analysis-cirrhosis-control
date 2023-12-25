# scRNA-seq-analysis-cirrhosis-control

Cirrhosis is a chronic liver disease that results from various etiologies, including viral
infections, alcohol abuse, and nonalcoholic steatohepatitis. The disease is characterized by the formation of fibrous scar tissue within the liver. This fibrosis is a significant contributor to liver damage and complications associated with the disease,
ultimately leading to Hepatocellular carcinoma (HCC). Single cell RNA sequencing
provides a valuable tool to study the gene expression patterns at a single cell level.
In this study, I aimed to compare the gene expression levels of established cirrhosis
and fibrosis markers between cirrhotic and normal cell samples using publicly available scRNA seq data from the GEO database. I have analyzed the obtained data
using the Seurat R package. After quality control, normalization, and other standard
analysis steps. The cluster cell types were annotated both manually using the
PanglaoDB, which includes data from various species, tissues, and cell types and
by using SingleR package, which provides a set of built-in references. The discrepancies in the annotations are viewed.


The single cell transcriptomic data was obtained from GEO
database with ascension number: GSE136103 (paper: Lev-
eraging SC-RNA-Seq DATA to uncover the association
between cell type and chronic liver disease). Three healthy
samples and three cirrhotic samples were used. The tran-
scripts were sequenced using Illumina Hiseq 4000 and all
the cirrhotic samples were collected from female patients
with the cause of disease being “NAFLD”.


Seurat Package (version 4.3.0) was used for data analysis
while SingleR package (Version 2.0.0) was used for cell
type annotation and Harmony package (Version 0.1.1) was
used for sample data integration

library(Seurat)
library(dplyr)
ovary_count<-Read10X("/media/disk1/RNA-seq/ovary_aggr/outs/filtered_feature_bc_matrix/")
ovary<- CreateSeuratObject(counts = ovary_count, project = "ovary", min.cells = 3, min.features = 200)

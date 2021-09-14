***************QC and clustering************************
library(Seurat)
library(dplyr)
library(harmony)
ovary_count<-Read10X("/media/disk1/RNA-seq/ovary_aggr/outs/filtered_feature_bc_matrix/")
ovary<- CreateSeuratObject(counts = ovary_count, project = "ovary", min.cells = 3, min.features = 200)
ovary[["percent.mt"]] <- PercentageFeatureSet(ovary, pattern = "^MT-")
ovary<- subset(ovary, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
ovary<- NormalizeData(ovary, normalization.method = "LogNormalize", scale.factor = 10000)
ovary<- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2100)
all.genes <- rownames(ovary)
ovary<- ScaleData(ovary, features = all.genes)
ovary<- RunPCA(ovary, features = VariableFeatures(object = ovary))


***************clustering************************
library(Seurat)
library(dplyr)
library(harmony)
ovary_count<-read.table("ovary_count",header=T,row.name=1,sep="\t")
ovary<- CreateSeuratObject(counts = ovary_count, project = "ovary")
ovary[["percent.mt"]] <- PercentageFeatureSet(ovary, pattern = "^MT-")
ovary<- NormalizeData(ovary, normalization.method = "LogNormalize", scale.factor = 10000)
ovary<- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2100)
all.genes <- rownames(ovary)
ovary<- ScaleData(ovary, features = all.genes,vars.to.regress = "percent.mt")
ovary<- RunPCA(ovary, features = VariableFeatures(object = ovary))
ovary@metadata$age<-ovary_metadata$age
ovary@metadata$group<-ovary_metadata$group
ovary <- ovary %>% RunHarmony("age", plot_convergence = TRUE)
ovary <- ovary %>% RunUMAP(reduction = "harmony", dims = 1:15) %>% FindNeighbors(reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.1) %>% identity()
DimPlot(ovary,pt.size = .5)
DimPlot(ovary,pt.size = .5,split.by="group")
DimPlot(ovary,pt.size = .5,split.by="age")
ovary_roc.markers <- FindAllMarkers(ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "roc")
new.cluster.ids <- c("Stromal_cell", "Blood_endothelial_cell", "Granulosa_cell", "Smooth_muscle_cell", "Immune_cell", "Lymphatic_endothelial_cell", "Epithelial_cell", "Theca_cell","Stromal_cell")
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5,split.by="age")
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5,split.by="group")
DotPlot(ovary,features = c("DCN","COL6A3","LUM","PDGFRA","VWF","FLT1","CDH2","SERPINE2","CYP19A1","INHA","FOXL2","AMH","MYH11","ACTA2","TAGLN","MCAM","PTPRC","CD53","CXCR4","PROX1","FLT4","PAX8","CDH1","CLDN1","STAR","CYP17A1"),cols = c("white","red"))
SC_aging_markers <- FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Stromal_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
BEC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Blood_endothelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
GC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Granulosa_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
SMC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Smooth_muscle_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
IC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Immune_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
LEC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Lymphatic_endothelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
EpiC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Epithelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
TC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Theca_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
write.table(SC_aging_markes,"SC_aging_markers",sep="\t")
write.table(SC_aging_markes,"BEC_aging_markers",sep="\t")
write.table(SC_aging_markes,"GC_aging_markers",sep="\t")
write.table(SC_aging_markes,"SMC_aging_markers",sep="\t")
write.table(SC_aging_markes,"IC_aging_markers",sep="\t")
write.table(SC_aging_markes,"LEC_aging_markers",sep="\t")
write.table(SC_aging_markes,"EpiC_aging_markers",sep="\t")
write.table(SC_aging_markes,"TC_aging_markers",sep="\t")

SC_aging_logfc<-FindMarkers(ovary,subset.ident = "Stromal_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
BEC_aging_logfc<-FindMarkers(ovary,subset.ident = "Blood_endothelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
GC_aging_logfc<-FindMarkers(ovary,subset.ident = "Granulosa_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
SMC_aging_logfc<-FindMarkers(ovary,subset.ident = "Smooth_muscle_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
IC_aging_logfc<-FindMarkers(ovary,subset.ident = "Immune_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
LEC_aging_logfc<-FindMarkers(ovary,subset.ident = "Lymphatic_endothelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
EpiC_aging_logfc<-FindMarkers(ovary,subset.ident = "Epithelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
TC_aging_logfc<-FindMarkers(ovary,subset.ident = "Theca_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
write.table(SC_aging_logfc,"SC_aging_logfc",sep = "\t")
write.table(BEC_aging_logfc,"BEC_aging_logfc",sep = "\t")
write.table(GC_aging_logfc,"GC_aging_logfc",sep = "\t")
write.table(SMC_aging_logfc,"SMC_aging_logfc",sep = "\t")
write.table(IC_aging_logfc,"IC_aging_logfc",sep = "\t")
write.table(LEC_aging_logfc,"LEC_aging_logfc",sep = "\t")
write.table(EpiC_aging_logfc,"EpiC_aging_logfc",sep = "\t")
write.table(TC_aging_logfc,"TC_aging_logfc",sep = "\t")
SC_aging_logfc<-SC_aging_logfc[,"avg_log2FC",drop=F]
BEC_aging_logfc<-BEC_aging_logfc[,"avg_log2FC",drop=F]
GC_aging_logfc<-GC_aging_logfc[,"avg_log2FC",drop=F]
SMC_aging_logfc<-SMC_aging_logfc[,"avg_log2FC",drop=F]
IC_aging_logfc<-IC_aging_logfc[,"avg_log2FC",drop=F]
LEC_aging_logfc<-LEC_aging_logfc[,"avg_log2FC",drop=F]
EpiC_aging_logfc<-EpiC_aging_logfc[,"avg_log2FC",drop=F]
TC_aging_logfc<-TC_aging_logfc[,"avg_log2FC",drop=F]
ovary_aging_logfc<-merge(SC_aging_logfc,BEC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,GC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,SMC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,IC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,LEC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,EpiC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,TC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
colnames(ovary_aging_logfc)<-c("SC","BEC","GC","SMC","IC","LEC","EpiC","TC")
write.table(ovary_aging_logfc,"ovary_aging_logfc",sep = "\t")

***************ovary_aging_DEG_heatmap************************
library(pheatmap)
library(RColorBrewer)
DEG_list<-read.table("DEG_list.txt",sep = "\t",row.names = 1,header = T)
ovary_aging_DEG<-merge(ovary_aging_logfc,DEG_list,by=0)
rownames(ovary_aging_DEG)<-ovary_aging_DEG$Row.names
ovary_aging_DEG<-ovary_aging_DEG[,2:9]
ovary_aging_DEG<-ovary_aging_DEG[rownames(DEG_list),]
pheatmap(as.matrix(ovary_aging_DEG),scale = "none",cluster_rows = F,cluster_cols = F,annotation_row = DEG_list,filename = "ovary_aging_DEG_heatmap.pdf",color=colorRampPalette(c("blue", "white", "red"))(200))

***************Common_DEG_heatmap************************
Common_up<-ovary_aging_logfc[grepl("^BCL6$|^SYNE1$|^SETBP1$|^SSBP2$|^RORA$|^CCNY$|^KMT2E$|^CHST11$|^EXOC4$|^KIF13A$|^IMMP2L$|^ALPK1$|^COG5$|^PLA2R1$|^ZCCHC7$|^ARID5B$|^XIST$|^NCOA1$|^ZBTB16$|^ANKRD36$|^CAMK2D$|^ZHX2$|^STAG1$|^KIAA1328$|^DYNC2H1$|^CUX1$|^OGA$|^NR3C2$|^AKT3$|^ZBTB20$|^PHF21A$|^MEIS2$|^MAST4$|^RALGAPA1$|^TNRC6B$|^NFIB$|^ZNF292$|^DENND4C$|^NCOA2$|^SENP6$|^CCND3$|^TNS1$|^MYCBP2$|^MLLT6$|^GULP1$|^CRADD$|^VPS13C$|^POGZ$|^USP53$|^SBF2$|^MIR29B2CHG$|^INPP4B$|^AGAP1$|^PTGIS$|^CDKAL1$|^HSPA1A$|^ANKS1B$|^L3MBTL4$|^VWA8$|^SLC9A9$|^PPP6R3$|^BABAM2$|^IGF1R$|^RNF149$|^ATP8A1$|^CHPT1$|^NAMPT$|^TMOD1$|^BCAS3$|^RUNX1T1$|^MYO1D$|^PITPNC1$|^PATJ$|^SLC25A37$|^RABGAP1L$|^NFIA$|^RNF217$|^ARID4B$|^DPYD$|^SNX9$|^DENND4A$|^ABLIM1$|^MKLN1$|^FARS2$|^DDX17$|^AKAP13$|^NBEA$|^USP24$|^WDPCP$|^KANSL1$|^C1orf21$|^LRMDA$|^GPCPD1$|^TCF7L2$|^DIAPH2$|^ATG7$|^ANAPC16$|^TACC1$|^SND1$|^ZEB1$|^SPPL3$|^PRRX1$|^FAF1$|^MLXIP$|^TBL1X$|^LRBA$|^PRKAG2$|^MAP3K5$|^FAM13A$|^HIPK3$|^TRAF3IP2-AS1$|^PPP3CC$|^SCAPER$|^NUBPL$|^TRAPPC9$|^NPAS3$|^RBMS3$|^BBS9$|^PCM1$|^STX8$|^VTI1A$|^ACSS3$|^TTC37$|^SPIDR$|^ADD3$|^LYST$|^VGLL4$|^FGD4$|^CBLB$|^PARD3$|^FNBP1$|^ARID1B$|^GPHN$|^BTAF1$|^SOS1$|^ABI3BP$|^DANT2$|^SPATA6$|^RBPMS$|^PLCG2$|^TMTC1$|^LIMK2$|^PRKRIP1$|^CLSTN2$|^ANKRD12$|^MAGI2$|^FBXO32$|^ELL2$|^ECHDC2$|^MAP4K5$|^N4BP2L1$|^GPATCH8$|^C3$|^CCSER1$|^C2CD5$|^UTRN$|^TPCN1$|^PTPRM$|^UMAD1$|^FTO$|^TBC1D22A$|^PDSS2$|^BOC$|^SOX5$|^NPEPPS$|^RICTOR$|^PLSCR1$|^SVIL$|^LARGE1$|^OXR1$|^NEK7$|^STK3$|^TMEM132C$|^C9orf3$|^KDM3A$|^JAZF1$|^LAMA2$|^RALGAPA2$|^EPHA6$|^NDRG1$|^CFLAR$|^KDM4C$|^GABARAPL1$|^EVA1C$|^FOXP1$|^C1RL$|^ELMO1$|^PDK1$|^ARHGAP29$|^TNXB$|^AFF1$|^MOB3B$|^ELP2$|^MAP3K20$|^MICAL3$|^ARHGAP26$|^MLLT10$|^TNRC6A$|^EFEMP1$|^AUH$|^GLIS3$|^CLASP1$|^NHS$|^UGP2$|^SERINC5$|^PRELID2$|^FCHO2$|^THADA$|^CDK19$|^CHST15$|^PGM5$|^FOXN3$|^CLASP2$|^FYN$|^BMPR2$|^RBMS1$|^MED13L$|^C9orf72$",rownames(ovary_aging_logfc)),]
Common_down<-ovary_aging_logfc[grepl("^IGFBP7$|^SPARC$|^TMSB4X$|^LINC00486$|^APOE$|^B2M$|^S100A6$|^TAGLN$|^JUND$|^JUN$|^TMSB10$|^RBP1$|^MTRNR2L12$|^MALAT1$|^RPL15$|^DHFR$|^MT-CO3$|^EIF1$|^SERF2$|^SRP14$|^INO80D$|^FDX1$|^DUSP1$|^VIM$|^MT-ATP6$|^RPL10A$|^RPL13A$|^RPS23$|^RPS19$|^BEX3$|^RPLP1$|^MT-CO2$|^RPL30$|^FAU$|^RPS2$|^RPL41$|^PEBP1$|^PSAP$|^ACTB$|^MYL12A$|^RPL11$|^APOC1$|^STAR$|^LINGO1$|^RASGEF1B$|^RPL3$|^RPS7$|^RPS14$|^MT-ND5$|^RPL13$|^CD63$|^RPL27$|^RPS8$|^FN1$|^RPL10$|^HBB$|^CYP11A1$|^KTN1$|^RPL35A$|^RPL7A$|^MT-ND1$|^COL1A2$|^RPS15$|^RPS18$|^MT-CO1$|^RPL18$|^MYL6$|^GSTA1$|^RPS29$|^CEBPB$|^LGALS1$|^RPS13$|^SESN3$|^MT-ND4L$|^RPS9$|^RPL28$|^KCNIP4$|^MT-ND2$|^RDX$|^BAIAP2L1$|^P3H2$|^MT-ND3$|^ANXA5$|^CST3$|^COL4A2$|^RAPH1$|^ENAH$|^EPS8$|^RPS4X$|^TPT1$|^SLC38A2$|^TCEAL4$|^RPL36$|^HSPD1$|^MT-CYB$|^OSBPL1A$|^FTL$|^CHCHD2$|^PRSS23$|^RPL29$|^PDXDC1$|^GAPDH$|^PSD3$|^TXNRD1$|^ACTG1$|^SORBS2$|^RPL14$|^RPS5$|^TCEAL9$|^FILIP1L$|^MT-ND4$|^MCU$|^ATP5F1E$|^GRAMD1B$|^SOX4$|^RPL27A$|^RPS3$|^COL1A1$|^ITPR1$|^PDE7B$|^RPL21$|^PTMS$|^S100A16$|^SF1$|^PPIA$|^TIMP1$|^ARHGAP18$|^ANGPT2$|^SCARB1$|^DOCK10$|^ACSM1$|^APOO$|^COX4I1$|^RPS25$|^RPL23A$|^RPS27A$|^ACSM3$|^MYL12B$|^GREB1$|^CD81$|^XPO5$|^RPS6$|^UBA52$|^AMBRA1$|^LAPTM4A$|^PPARG$|^RPL19$|^NACA$|^RPL5$|^HIBCH$|^COL4A1$|^NEAT1$|^DDX24$|^PRKCA$|^CPEB4$|^MTRNR2L8$|^MEF2C$|^RPS24$|^RAB11A$|^HSP90AB1$|^UBC$|^SMARCA1$|^GATA4$|^S100A13$|^OOEP$|^HMGB1$|^RPL23$|^HSPA4$|^RPL6$|^C11orf96$|^CBR4$|^PTMA$|^LRP4$|^ATRNL1$|^THUMPD1$|^EXT1$|^MIB1$|^GJA1$|^RPL7$|^MYL9$|^GPC6$|^RACK1$",rownames(ovary_aging_logfc)),]
Common_DEG<-rbind(Common_up,Common_down)
pheatmap(as.matrix(Common_DEG),scale = "none",cluster_rows = T,cluster_cols = T,filename = "Common_DEG_heatmap.pdf",color=colorRampPalette(c("blue", "white", "red"))(200))

***************cell_specific_DEG_heatmap************************
cell_specific_DEG<-read.table("cell_specific_DEG",sep = "\t",row.names = 1,header = T)
ovary_cell_specific_DEG<-merge(ovary_aging_logfc,cell_specific_DEG,by=0)
rownames(ovary_aging_DEG)<-ovary_aging_DEG$Row.names
ovary_aging_DEG<-ovary_aging_DEG[,2:9]
ovary_aging_DEG<-ovary_aging_DEG[rownames(DEG_list),]
pheatmap(as.matrix(ovary_aging_DEG),scale = "none",cluster_rows = F,cluster_cols = F,annotation_row = DEG_list,filename = "ovary_aging_DEG_heatmap.pdf",color=colorRampPalette(c("blue", "white", "red"))(200))

#2022Sep27, Ting Sun
#DPZ Winkler data
#cell culture bulk RNA-seq
#deconvolution analysis

#mostly used MuSiC2 pipeline
#for the purpose of studying the treated cells
#ref:https://xuranw.github.io/MuSiC/articles/pages/MuSiC2.html

#run under r4-base

library(MuSiC)

library(biomaRt)
library(Biobase)

library(dplyr)
library(tidyr)

#first assemble the bulk expression set
exp<-read.csv("indir/DPZ_Winkler_bulkRNA_allsamples_RAW_value.csv", 
              stringsAsFactors = FALSE)

head(exp)

rownames(exp)<-exp$Gene_symbol
exp<-exp[,-c(1:3)]

head(exp)
dim(exp)

#transfer genes to ensemble ID
#to fit scRNA-seq data for deconvolution
anno<-read.csv("refdir/Gene_translations/Mmul10_to_hsapiens_annotation_ENSEMBLE_GeneSymbol.csv")

head(anno)

exp$ENS<-NA

for(i in 1:nrow(exp)){
    position<-match(as.character(rownames(exp))[i], anno$Gene.name)[1]
    exp$ENS[i]<-as.character(anno$Gene.stable.ID.1)[position]
}

head(exp)

#still need to rerun alignment
#at the moment proceed with filter out nNA entries

exp<-na.omit(exp)

length(unique(exp$ENS))

dim(exp)

exp_unique<- exp[ !duplicated(exp$ENS), ]

rownames(exp_unique)<-exp_unique$ENS

head(exp_unique)

exp_unique<-exp_unique[,-3]

head(exp_unique)

meta<-data.frame(row.names = c("Sample1_raw","Sample2_raw"))
meta$Species<-"Macaques"
meta$Tissue<-"MixCulture"
meta$Groups<-rownames(meta)

meta

phenoData =  new("AnnotatedDataFrame", data = meta)

bulk.set<-ExpressionSet(assayData = data.matrix(exp_unique), phenoData =  phenoData)

bulk.set

#saveRDS(bulk.set,
#       "indir/2022DPZ_Winkler_bulkRNA_HumanGeneTranslated_ExpSet.rds")

########################################
#script break
#following are deconvolution pipeline
##########################################

bulk.set<-readRDS("indir/2022DPZ_Winkler_bulkRNA_HumanGeneTranslated_ExpSet.rds")

bulk.set

#read in deconvolution reference
sc.sce<-readRDS("indir/external_RAW/scRNA/GSE127774_MacaquesACC_scExp.rds")

sc.sce

#confirm deconvolution groups
table(bulk.set$Groups)

#fix x array bug during deconcolution
counts(sc.sce) <- as.matrix(counts(sc.sce))
logcounts(sc.sce) <- as.matrix(logcounts(sc.sce))

bulk.mtx = exprs(bulk.set)

est.prop<-music_prop(bulk.mtx = bulk.mtx, sc.sce = sc.sce,
                     clusters = "CellType", samples = "SampleID")

names(est.prop)

est.prop

#saveRDS(est.prop,
#       "Outdir/2022_DPZ_Winkler_RNAseq_deconvolution_result.rds")

########################################################################

###############script break
#from here deconvolution result is already calculated and saved
#recall object to run the following script

est.prop<-readRDS("Outdir/2022_DPZ_Winkler_RNAseq_deconvolution_result.rds")

#read in annotation file
anno<-read.csv("refdir/Gene_translations/Mmul10_to_hsapiens_annotation_ENSEMBLE_GeneSymbol.csv")

sum(est.prop$Est.prop.weighted[1,])

#test visualizations
#stack barplot
sample_per_genotype<-as.data.frame(est.prop$Est.prop.weighted)
    sample_per_genotype

sample_per_genotype$sampleID<-rownames(sample_per_genotype)

sample_per_genotype

#prepare plotting data frame
(df<-pivot_longer(sample_per_genotype, cols =1:6, names_to = "Celltype", values_to = "Freq"))

(plot1<-ggplot(df, aes(y=Freq, x= sampleID, 
                                  fill= Celltype)) +  
                geom_bar( stat="identity", position="fill", width=0.6)+
                theme(axis.text.x = element_text(size = 20, #angle = 45, 
                                 hjust = 0.4, vjust = 0.5),
                  axis.text.y = element_text(size = 15),
                     axis.title.x = element_blank(),
                    legend.text = element_text(size=15),
                     legend.title = element_text(size=15))+
                ylab("")+
                scale_x_discrete(labels= unique(sample_per_genotype$sampleID))+ 
                coord_cartesian(expand = FALSE)+
                scale_fill_manual(values=c("lightblue3","gray90","rosybrown2",
                                           "tan1","#91D1C2B2","lightcoral"))+
                theme(legend.position="right",
                     panel.background = element_blank()))

















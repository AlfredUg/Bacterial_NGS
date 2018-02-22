
##############Classification results ##########################
abund_table <- read.csv("~/Dropbox/Msc_Research_Project/Bacterial_NGS/final_results/abund_table.csv",row.names = 1)
333333333333333333333333
grouping_info <- colnames(abund_table)

bac_samples <- strsplit(grouping_info, "_")
bac_samples
sampleType<-NULL
for(i in 1:length(bac_samples)){
  sampleType[i]<-bac_samples[[i]][2]
}

require(phyloseq)
SAM<-sample_data(data.frame(row.names = colnames(abund_table), Location=as.character(sampleType)))
SAM<-SAM[which(!(SAM$Location%in%c("S8", "S9", "S10"))),]
abund_table <- t(abund_table)
abund_table <- abund_table[which(rownames(abund_table)%in%rownames(SAM)),]
OTU<-otu_table(abund_table, taxa_are_rows = F )
phys<-phyloseq(OTU, SAM)

require(microbiomeSeq)
require(gtable)
require(gridExtra)
require(ggplot2)
require(grid)

p1 <- plot_anova_diversity(phys, method = c("richness","simpson","shannon"), grouping_column = "Location", pValueCutoff = 0.05)
p1 <- p1 + guides(fill=guide_legend(title="Location")) + xlab("Location")
p1

physp<-normalise_data(phys, norm.method = "relative")
pp<-plot_taxa(physp, grouping_column = "Location", filename = "AFI.csv")
pp <- pp+ylab("Relative abundance")

ord.res <- ordination(phys, method="PCoA", grouping_column="Location")
p2<-plot.ordination(ord.res, method="PCoA", pvalue.cutoff = 0.05, show.pvalues = T, N = 5, extra_marginspace = 0.08)
p2


ord.res1 <- ordination(phys, method="NMDS", grouping_column="Location")
p3<-plot.ordination(ord.res1, method="NMDS", pvalue.cutoff = 0.05, show.pvalues = T, N = 5, extra_marginspace = 0.35)
p3

deseq_sig <- differential_abundance(phys, grouping_column = "Location",output_norm = "log-relative")
p4<-plot_signif(deseq_sig$plotdata, top.taxa = 15)
p4 <- p4 + ylab("Log-relative normalised abundance") + scale_colour_manual("Location", values = c("red","green","blue")) + xlab("Location")
p4

phys<-normalise_data(phys, norm.method = "log-relative")
krus_sig <- kruskal_abundance(phys, grouping_column = "Location", pvalue.threshold = 0.05)

####################### plots for virulence factor results #################
above90vfs <- read.delim("~/Dropbox/Msc_Research_Project/Bacterial_NGS/final_results/virulence/above90vfs.txt", header=FALSE)
above90vfs <- above90vfs[,-c(2,4)]
colnames(above90vfs) <- c("AccIds","VirulenceFactor","Pathogen" ,"PercentMatch", "Evalue", "Sample")

# 1. GET NUMBER OF VIRULENCE FACTORS FOR EACH PATHOGEN
pathogens <- unique(above90vfs$Pathogen)
df<-NULL
for(i in pathogens){
  nv<-length(unique(above90vfs$VirulenceFactor[above90vfs$Pathogen==i]))
  tmp<-data.frame("Pathogen"=i, "No_VFs"=nv)
  if(is.null(df)){
    df<-tmp
  }
  else{
    df <- rbind(df, tmp)
  }
}
#print(df)
p <- ggplot(data = df,aes(x=Pathogen,y=No_VFs)) + theme_bw()
p <- p+geom_bar(stat = "identity",position=position_dodge(), width = 0.5, fill="darkblue")
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p <- p + ylab("Number of virulence factors")
p
# 2. GET THE DISTRIBUTION OF VFs IN THE DIFFERENT SAMPLE LOCATIONS
vfs_samples <- strsplit(gsub("-","_", as.character(above90vfs$Sample)), "_")
sampleType<-NULL
for(i in 1:length(vfs_samples)){
  sampleType[i]<-vfs_samples[[i]][2]
}
above90vfs$SampleType <- sampleType

df1<-NULL
for(i in unique(above90vfs$SampleType)){
  samplePathogens<-unique(above90vfs$Pathogen[above90vfs$SampleType==i])
  for(j in samplePathogens){
    nps<- length(unique(above90vfs$Pathogen)[above90vfs$Pathogen==j & above90vfs$SampleType==i])
    tmp<-data.frame("Pathogen"=j, "No_Pathogens"=nps, "SampleType"=i)
    if(is.null(df1)){
      df1<-tmp
    }
    else{
      df1 <- rbind(df1, tmp)
    }
  }
}
#print(df1)
p <- ggplot(data = df1,aes(x=Pathogen,y=No_Pathogens, fill=SampleType)) + theme_bw()
p <- p+geom_bar(stat = "identity",position=position_dodge(), width = 0.5)
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ylab("Number of Virulence factors") + xlab("Potential pathogens")
p <- p + guides(fill=guide_legend(title="Sample Location"))
p

vf_pathogens<-NULL

for(i in unique(above90vfs$Pathogen)){
  no_class_1 <- length(unique(above90vfs$Pathogen)[above90vfs$Pathogen==i & above90vfs$SampleType=="1"])
  no_class_2 <- length(unique(above90vfs$Pathogen)[above90vfs$Pathogen==i & above90vfs$SampleType=="2"])
  no_class_3 <- length(unique(above90vfs$Pathogen)[above90vfs$Pathogen==i & above90vfs$SampleType=="3"])
  
  tmp<-data.frame("Class"=i, "1"=no_class_1, "2"=no_class_2, "3"=no_class_3)
  
  if(is.null(vf_pathogens)){
    vf_pathogens<-tmp
  }
  else{
    vf_pathogens <- rbind(vf_pathogens, tmp)
  }
}
vf_pathogens

ch<-above90vfs[,c("VirulenceFactor", "Pathogen")]
ch_rdp<- ch[!(duplicated(ch)),]

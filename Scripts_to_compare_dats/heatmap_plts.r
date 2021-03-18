# make heat maps for paper
# using OTU tables from Batya

library(ggplot2)
library(reshape2)
library(RColorBrewer)
setwd('./Desktop/Old projects/everestWGS/OTUtables_makeheatmaps/')

# read in files
OTUtables <- list.files(path='.', pattern='.csv', full.names=T)
Data1_OTU_wgs_16s <- read.csv(OTUtables[1], header=T)
Data2_OTU_wgs_microb <- read.csv(OTUtables[2], header=T)
Data3_OTU_metabar <- read.csv(OTUtables[3], header=T)
Data5_OTU_wgs_refmapped <- read.csv(OTUtables[4], header=T)

# format tables - from matrix to df
data1 <- melt(Data1_OTU_wgs_16s); colnames(data1) <- c('Site', 'Order', 'NumReads')
data2 <- melt(Data2_OTU_wgs_microb); colnames(data2) <- c('Site', 'Order', 'NumReads')
data3 <- melt(Data3_OTU_metabar); colnames(data3) <- c('Site', 'Order', 'NumReads')
data5 <- melt(Data5_OTU_wgs_refmapped); colnames(data5) <- c('Site', 'Order', 'NumReads')

# as bar plots
order.cols.dat1 <- length(unique(data1$Order))
ordercolors.dat1 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat1)
order.cols.dat2 <- length(unique(data2$Order))
ordercolors.dat2 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat2)
order.cols.dat3 <- length(unique(data3$Order))
ordercolors.dat3 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat3)
order.cols.dat5 <- length(unique(data5$Order))
ordercolors.dat5 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat5)

ggplot(data1, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat1) +
  ggtitle('Data1_OTU_wgs_16s') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('Data1_OTU_wgs_16s_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data2, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat2) +
  ggtitle('Data2_OTU_wgs_microb') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('Data2_OTU_wgs_microb_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data3, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat3) +
  ggtitle('Data3_OTU_metabar') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('Data3_OTU_metabar_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data5, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat5) +
  ggtitle('Data5_OTU_wgs_refmapped') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('Data5_OTU_wgs_refmapped_barplot.jpg', height=8, width=20, units='in', dpi=600)

# as heat maps
ggplot(data1, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='black') +
  scale_fill_viridis_c(direction=-1,
                       limits=c(0, max(data1$NumReads))) +
  ggtitle('Data1_OTU_wgs_16s') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave('Data1_OTU_wgs_16s_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data2, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='black') +
  scale_fill_viridis_c(direction=-1,
                       limits=c(0, max(data2$NumReads))) +
  ggtitle('Data2_OTU_wgs_microb') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave('Data2_OTU_wgs_microb_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data3, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='black') +
  scale_fill_viridis_c(direction=-1,
                       limits=c(0, max(data3$NumReads))) +
  ggtitle('Data3_OTU_metabar') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave('Data3_OTU_metabar_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data5, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='black') +
  scale_fill_viridis_c(direction=-1,
                       limits=c(0, max(data5$NumReads))) +
  ggtitle('Data5_OTU_wgs_refmapped') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave('Data5_OTU_wgs_refmapped_heatmap.jpg', height=15, width=10, units='in', dpi=600)

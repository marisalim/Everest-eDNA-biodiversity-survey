# make heat maps for paper
# using OTU tables from Batya

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
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

# clean up files
# data1 should not have Stamenopiles (euk) for the singlem results
data1_clean <- data1[data1$Order!='Stramenopiles', ]

# ---- as bar plots ----
order.cols.dat1 <- length(unique(data1_clean$Order))
ordercolors.dat1 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat1)
order.cols.dat2 <- length(unique(data2$Order))
ordercolors.dat2 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat2)
order.cols.dat3 <- length(unique(data3$Order))
ordercolors.dat3 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat3)
order.cols.dat5 <- length(unique(data5$Order))
ordercolors.dat5 <- colorRampPalette(brewer.pal(8, "Set1"))(order.cols.dat5)

ggplot(data1_clean, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat1) +
  scale_y_continuous(labels=scales::comma) +
  ggtitle('Data1_OTU_wgs_16s') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data1_OTU_wgs_16s_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data2, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat2) +
  scale_y_continuous(labels=scales::comma) +
  ggtitle('Data2_OTU_wgs_microb') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data2_OTU_wgs_microb_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data3, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat3) +
  scale_y_continuous(labels=scales::comma) +
  ggtitle('Data3_OTU_metabar') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data3_OTU_metabar_barplot.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data5, aes(x=Site, y=NumReads, fill=factor(Order))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=ordercolors.dat5) +
  scale_y_continuous(labels=scales::comma) +
  ggtitle('Data5_OTU_wgs_refmapped') +
  ylab('Number of seq reads') +
  xlab('Sample') +
  guides(fill=guide_legend(title="Order")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data5_OTU_wgs_refmapped_barplot.jpg', height=8, width=20, units='in', dpi=600)

# ---- as heat maps ----
# make 0 counts NA's so they can be plotted as white boxes in tile plots
data1_clean[data1_clean == 0] <- NA
data2[data2==0] <- NA
data3[data3==0] <- NA
data5[data5==0] <- NA

ggplot(data1_clean, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='grey') +
  scale_fill_gradientn(colors=viridis_pal()(20),
                       limits=c(0, max(data1_clean$NumReads)),
                       na.value='white',
                       labels=comma) +
  ggtitle('Data1_OTU_wgs_16s') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data1_OTU_wgs_16s_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data2, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='grey') +
  scale_fill_gradientn(colors=viridis_pal()(20),
                       limits=c(0, max(data2$NumReads)),
                       na.value='white',
                       labels=comma) +
  ggtitle('Data2_OTU_wgs_microb') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data2_OTU_wgs_microb_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data3, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='grey') +
  scale_fill_gradientn(colors=viridis_pal()(20),
                       limits=c(0, max(data3$NumReads)),
                       na.value='white',
                       labels=comma) +
  ggtitle('Data3_OTU_metabar') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data3_OTU_metabar_heatmap.jpg', height=8, width=10, units='in', dpi=600)

ggplot(data5, aes(x=Site, y=factor(Order))) +
  geom_tile(aes(fill=NumReads), col='grey') +
  scale_fill_gradientn(colors=viridis_pal()(20),
                       limits=c(0, max(data5$NumReads)),
                       na.value='white',
                       labels=comma) +
  ggtitle('Data5_OTU_wgs_refmapped') +
  ylab('Taxonomic Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=15, face='bold'))
ggsave('Data5_OTU_wgs_refmapped_heatmap.jpg', height=15, width=10, units='in', dpi=600)

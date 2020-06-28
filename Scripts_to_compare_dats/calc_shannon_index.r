# shannon's index and simpson's index

library(vegan)
library(ggplot2)
library(cowplot)
library(reshape2)

# ---- input WGS data ----
## this version is just the Kraken relative abundance %
genewiz_relabund <- read.csv('../../Volumes/Seagate Backup Plus Drive/WCS_MarisaLimpostdocBACKUP/ALL_eDNA_project_materials/Everest_metagen_dat/30-254524058_genewiz_metagenomics_analysis/projectData/data/30-254524058_nozeros_otu.csv')
genewiz_taxaIDs <- read.csv('../../Volumes/Seagate Backup Plus Drive/WCS_MarisaLimpostdocBACKUP/ALL_eDNA_project_materials/Everest_metagen_dat/30-254524058_genewiz_metagenomics_analysis/projectData/data/30-254524058_taxonomy.csv', header=TRUE)
colnames(genewiz_taxaIDs) <- c('TaxID', 'Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
newdat <- merge(genewiz_relabund, genewiz_taxaIDs, by='TaxID')

# ---- WGS OTU counts per taxa rank ----
# hardly any rows if subset to only Domain or Kingdom rank, so these are not useful

# phylum rank subset
toPhylum <- newdat[newdat$Phylum != '<NA>' & newdat$Class == '<NA>' &  newdat$Order == '<NA>' &
                     newdat$Family == '<NA>' & newdat$Genus == '<NA>' & newdat$Species == '<NA>',]
# class rank subset
toClass <- newdat[newdat$Class != '<NA>' &  newdat$Order == '<NA>' &
                    newdat$Family == '<NA>' & newdat$Genus == '<NA>' & newdat$Species == '<NA>',]
# order rank subset
toOrder <- newdat[newdat$Order != '<NA>' & newdat$Family == '<NA>' & newdat$Genus == '<NA>' & newdat$Species == '<NA>',]
# family rank subset
toFamily <- newdat[newdat$Family != '<NA>' & newdat$Genus == '<NA>' & newdat$Species == '<NA>',]

OTUcountformatting <- function(OTUtab, therankcol){
  newdat_todf <- OTUtab[,c(therankcol,3:9,11:16,18:20,22:25)]
  rownames(newdat_todf) <- newdat_todf[,1]
  newdat_todf2<- newdat_todf[,c(-1)]
  newdat_df <- data.frame('OTUrankcounts'=rowSums(newdat_todf2))
  newdat_df$taxID <- row.names(newdat_df)
  return(newdat_df)
}
phy_df <- OTUcountformatting(toPhylum, 28)
class_df <- OTUcountformatting(toClass, 29)
order_df <- OTUcountformatting(toOrder, 30)
fam_df <- OTUcountformatting(toFamily, 31)

dim(phy_df) #23 phyla
dim(class_df) #59 classes
dim(order_df) #129 orders
dim(fam_df) #213 families

ggplot(phy_df, aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Total relative proportion of \nreads per taxon ID') +
  ggtitle('WGS OTU counts\nOTU = Phylum (n=23)') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUphylum_OTUcounts_genewizWGSredo.jpg',
       height=8, width=8, dpi=500)

ggplot(class_df, aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Total relative proportion of \nreads per taxon ID') +
  ggtitle('WGS OTU counts\nOTU = Class (n=59)') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUclass_OTUcounts_genewizWGSredo.jpg',
       height=8, width=10, dpi=500)

ggplot(order_df, aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Total relative proportion of \nreads per taxon ID') +
  ggtitle('WGS OTU counts\nOTU = Order (n=129)') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=5),
        axis.text.y=element_text(size=10))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUorder_OTUcounts_genewizWGSredo.jpg',
       height=8, width=15, dpi=500)

ggplot(fam_df, aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Total relative proportion of \nreads per taxon ID') +
  ggtitle('WGS OTU counts\nOTU = Family (n=213)') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=5),
        axis.text.y=element_text(size=10))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUfamily_OTUcounts_genewizWGSredo.jpg',
       height=8, width=20, dpi=500)

# ---- WGS filtered to Order rank ----
# now we need to do some filtering
# check if the remaining taxa in OTU table are the ones you want!
View(toOrder) # no more Family = NA rows
unique(toOrder$Domain) #Bacteria  Archaea   Eukaryota Viruses 

# save 'OTU' table
write.csv(toOrder[,c(30, 3:9,11:16,18:20,22:25)], 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGS_OTUbyOrder_table.csv', row.names=FALSE)

# calculate shannon's index
H_rough3 <- diversity(t(toOrder[,c(3:9,11:16,18:20,22:25)]), index='shannon')
dim(toOrder)
length(H_rough3)
D_rough3 <- diversity(t(toOrder[,c(3:9,11:16,18:20,22:25)]), index='invsimpson')

HD_rough_df3 <- data.frame('Shannon'=H_rough3, 'Simpson'=D_rough3)
HD_rough_df3$sampleID <- row.names(HD_rough_df3)

# add locality/site info to aggregate data by site
evsites <- rep('na', nrow(HD_rough_df3))
evsites[HD_rough_df3$sampleID=='TS13' |
          HD_rough_df3$sampleID=='TS14'] <- 'AboveEBC'
evsites[HD_rough_df3$sampleID=='TS15' |
          HD_rough_df3$sampleID=='TS17'] <- 'KalaPattarLake1'
evsites[HD_rough_df3$sampleID=='TS18' |
          HD_rough_df3$sampleID=='TS19'] <- 'KalaPattarLake2'
evsites[HD_rough_df3$sampleID=='TS20' |
          HD_rough_df3$sampleID=='TS21'] <- 'KongmaLaLake1_inlet'
evsites[HD_rough_df3$sampleID=='TS22' |
          HD_rough_df3$sampleID=='TS23'] <- 'KongmaLaLake1_outlet'
evsites[HD_rough_df3$sampleID=='TS24' |
          HD_rough_df3$sampleID=='TS25'] <- 'KongmaLaLake2'
evsites[HD_rough_df3$sampleID=='TS26' |
          HD_rough_df3$sampleID=='TS27'] <- 'KongmaLaLake3'
evsites[HD_rough_df3$sampleID=='TS28' |
          HD_rough_df3$sampleID=='TS29'] <- 'NuptseGlacierMtnStream'
evsites[HD_rough_df3$sampleID=='TS30' |
          HD_rough_df3$sampleID=='TS31'] <- 'KongmaLaLakeMtnStream'
evsites[HD_rough_df3$sampleID=='TS32' |
          HD_rough_df3$sampleID=='TS33'] <- 'LakeSouthNuptse'
HD_rough_df3$evsites <- evsites

HD_rough_df3_long <- melt(HD_rough_df3)

order_index <- ggplot(HD_rough_df3_long, aes(x=sampleID, y=value, fill=evsites)) +
  facet_wrap(~variable, ncol=1, scales='free') +
  geom_point(size=4, shape=21) +
  scale_fill_brewer(palette='Paired') +
  xlab('Sample ID') +
  ylab('Index value') +
  labs(fill='Sites') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        panel.grid.major=element_line(size=0.5, colour='seashell'),
        panel.grid.minor=element_line(size=0.5, colour='seashell'),
        legend.position='bottom',
        legend.text=element_text(size=8)) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))
order_otu <- ggplot(order_df[order(order_df$OTUrankcounts, decreasing = TRUE),][c(1:20),], aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Total relative proportion of \nreads per taxon ID') +
  # ggtitle('WGS OTU counts\nOTU = Order (top 20)') +
  ggtitle('Top 20 Order WGS OTUs') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10))
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUorder_shansimp_genewizWGSredo.jpg',
     height=8, width=15, units='in', res=500)
plot_grid(order_index, order_otu, ncol=2, align='v', labels=c('a)', 'b)'))
dev.off()

ggsave(plot=order_index, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/shansimp_genewizWGSredo.jpg',
       height=8, width=8, dpi=500)
# Counts by Order, per sample - grouped by site
#Order WGS OTUs with XX+ counts per sample
# give Order an ID number and we'll add a key to the name
otu_df <- data.frame(toOrder[,c(1, 30, 3:9,11:16,18:20,22:25)])
#add index number for the taxa (TaxID would be the actually taxa number, but not easy to
# read as a key for each bar, so let's try starting from 0)
dim(otu_df)
otu_df$index_for_fig <- c(1:nrow(otu_df))
wgsOrderdf <- melt(otu_df, id.vars=c('TaxID', 'Order', 'index_for_fig'))
count_thresh <- 0.2

wgsOrderdf$samporder <- factor(wgsOrderdf$variable,
                               levels=c('TS13', 'TS14', 'TS15', 'TS17',
                                        'TS18', 'TS19', 'TS20', 'TS21', 'TS22',
                                        'TS23', 'TS24', 'TS25', 'TS26', 'TS27',
                                        'TS28', 'TS29', 'TS30', 'TS31', 'TS32',
                                        'TS33'))

ggplot(wgsOrderdf, aes(x=Order, y=samporder, fill=value)) +
  geom_tile(alpha=0.5) +
  scale_fill_distiller(type='div', palette=9) +
  ylab('Sample ID') +
  xlab('Taxonomic Order') +
  labs(fill='Total relative proportion \nof reads per Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7), axis.text.y=element_text(size=8),
        legend.position='right', legend.text=element_text(size=8))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/OTUorderheatmap_genewizWGS.jpg',
     height=6, width=15, dpi=500)
     
# ----- MAIN FIGS: by site individually: USE THIS -----
write.csv(otu_df[,c(2,23,1,3:22)], 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGS_OrderOTUrelabund_key_data.csv', row.names=FALSE)
length(unique(wgsOrderdf[wgsOrderdf$value >=0.2,]$index_for_fig))# there are 29 Orders with >=0.2 counts
fig_orders <- unique(wgsOrderdf[wgsOrderdf$value >=0.2,c(2:3)])
fig_orders[order(fig_orders$index_for_fig),]

aboveEBC <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS13' | wgsOrderdf$variable=='TS14') & wgsOrderdf$value >= count_thresh, ], 
       aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Above Everest Base Camp \nTS13 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS13' & wgsOrderdf$value >= count_thresh, ]),
                ' TS14 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS14' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                             legend.position='right')
KPL1 <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS15' | wgsOrderdf$variable=='TS17') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kala Pattar Lake 1 \nTS15 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS15' & wgsOrderdf$value >= count_thresh, ]),
                ' TS17 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS17' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                        legend.position='right')
KPL2 <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS18' | wgsOrderdf$variable=='TS19') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kala Pattar Lake 2 \nTS18 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS18' & wgsOrderdf$value >= count_thresh, ]),
                ' TS19 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS19' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                        legend.position='right')
KLL1in <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS20' | wgsOrderdf$variable=='TS21') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 1 (inlet) \nTS20 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS20' & wgsOrderdf$value >= count_thresh, ]),
                ' TS21 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS21' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                                legend.position='right')
KLL1out <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS22' | wgsOrderdf$variable=='TS23') & wgsOrderdf$value >= count_thresh, ], 
                 aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 1 (outlet) \nTS22 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS22' & wgsOrderdf$value >= count_thresh, ]),
                ' TS23 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS23' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                                 legend.position='right')
KLL2 <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS24' | wgsOrderdf$variable=='TS25') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 2 \nTS24 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS24' & wgsOrderdf$value >= count_thresh, ]),
                ' TS25 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS25' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                      legend.position='right')
KLL3 <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS26' | wgsOrderdf$variable=='TS27') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 3 \nTS26 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS26' & wgsOrderdf$value >= count_thresh, ]),
                ' TS27 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS27' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                      legend.position='right')
NGMS <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS28' | wgsOrderdf$variable=='TS29') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Nuptse Glaciere Mtn Stream \nTS28 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS28' & wgsOrderdf$value >= count_thresh, ]),
                ' TS29 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS29' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                               legend.position='right')
KLLMS <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS30' | wgsOrderdf$variable=='TS31') & wgsOrderdf$value >= count_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake Mtn Stream \nTS30 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS30' & wgsOrderdf$value >= count_thresh, ]),
                ' TS31 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS31' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
        legend.position='right')
LSN <- ggplot(wgsOrderdf[(wgsOrderdf$variable=='TS32' | wgsOrderdf$variable=='TS33') & wgsOrderdf$value >= count_thresh, ], 
              aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Lake South Nuptse \nTS32 Orders n=', 
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS32' & wgsOrderdf$value >= count_thresh, ]),
                ' TS33 Orders n=',
                nrow(wgsOrderdf[wgsOrderdf$variable=='TS33' & wgsOrderdf$value >= count_thresh, ]),
                sep='')) + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                               legend.position='right')

jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrderOTUabove20_bysite1-5.jpg',
     height=20, width=15, units='in', res=500)
plot_grid(aboveEBC, KPL1, KPL2, KLL1in, KLL1out, ncol=1, align='h')
dev.off()
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrderOTUabove20_bysite6-10.jpg',
     height=20, width=15, units='in', res=500)
plot_grid(KLL2, KLL3, NGMS, KLLMS, LSN, ncol=1, align='h')
dev.off()

#try a different layout
panel1 <- plot_grid(aboveEBC, KPL1, KPL2, ncol=1)
panel2 <- plot_grid(KLL1in, KLL1out, KLL2, ncol=1)
panel3 <- plot_grid(KLL3, NGMS, KLLMS, ncol=1)
panel4 <- plot_grid(order_index, LSN, ncol=1, rel_heights = c(2,1))
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrder_shansimp_OTUabove20_bysite.jpg',
     height=15, width=25, units='in', res=500)
plot_grid(panel1, panel2, panel3, panel4, ncol=4)
dev.off()

# saving individual plots
ggsave(plot=aboveEBC, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_aboveEBC.jpg', height=6, width=10, dpi=500)
ggsave(plot=KPL1, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KPL1.jpg', height=6, width=10, dpi=500)
ggsave(plot=KPL2, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KPL2.jpg', height=6, width=10, dpi=500)
ggsave(plot=KLL1in, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KLL1in.jpg', height=6, width=10, dpi=500)
ggsave(plot=KLL1out, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KLL1out.jpg', height=6, width=10, dpi=500)
ggsave(plot=KLL2, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KLL2.jpg', height=6, width=10, dpi=500)
ggsave(plot=KLL3, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KLL3.jpg', height=6, width=10, dpi=500)
ggsave(plot=NGMS, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_NGMS.jpg', height=6, width=10, dpi=500)
ggsave(plot=KLLMS, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_KLLMS.jpg', height=6, width=10, dpi=500)
ggsave(plot=LSN, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/WGSOrders_LSN.jpg', height=6, width=10, dpi=500)

# ---- metabarcoding by Order ----
# this data input comes from running steps in the compare_eDNAdat_splits.r script
bysitedat <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/allEvdata_taxinfo_bysite.csv', header=TRUE)

# get just metabar contigs, also remove singletons
data_wide <- bysitedat[bysitedat$datatype=='metabar_contigs' &
                         bysitedat$conf > 1, ]
dim(bysitedat[bysitedat$datatype=='metabar_contigs',]) 
dim(data_wide) #none only have 1 read/contig, all were already 2+ reads/contig!

View(data_wide)
data_x <- table(data_wide$sampleID,
                data_wide$order)

# remove cols with only 0s
# no hits for ts32 metabarcoding, must remove too
datawide2 <- data_x[row.names(data_x) != 'ts16' & row.names(data_x) != 'ts35' & row.names(data_x) != 'ts32', colSums(data_x != 0) > 0]
dim(datawide2)

write.csv(t(datawide2), 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOTUbyOrder_table.csv')
H_met <- diversity(datawide2, index='shannon') 
D_met <- diversity(datawide2, index='invsimpson')

HD_met_df <- data.frame('Shannon'=H_met, 'Simpson'=D_met)
HD_met_df$sampleID <- row.names(HD_met_df)

evsites <- rep('na', nrow(HD_met_df))
evsites[HD_met_df$sampleID=='ts13' |
          HD_met_df$sampleID=='ts14'] <- 'AboveEBC'
evsites[HD_met_df$sampleID=='ts15' |
          HD_met_df$sampleID=='ts17'] <- 'KalaPattarLake1'
evsites[HD_met_df$sampleID=='ts18' |
          HD_met_df$sampleID=='ts19'] <- 'KalaPattarLake2'
evsites[HD_met_df$sampleID=='ts20' |
          HD_met_df$sampleID=='ts21'] <- 'KongmaLaLake1_inlet'
evsites[HD_met_df$sampleID=='ts22' |
          HD_met_df$sampleID=='ts23'] <- 'KongmaLaLake1_outlet'
evsites[HD_met_df$sampleID=='ts24' |
          HD_met_df$sampleID=='ts25'] <- 'KongmaLaLake2'
evsites[HD_met_df$sampleID=='ts26' |
          HD_met_df$sampleID=='ts27'] <- 'KongmaLaLake3'
evsites[HD_met_df$sampleID=='ts28' |
          HD_met_df$sampleID=='ts29'] <- 'NuptseGlacierMtnStream'
evsites[HD_met_df$sampleID=='ts30' |
          HD_met_df$sampleID=='ts31'] <- 'KongmaLaLakeMtnStream'
evsites[HD_met_df$sampleID=='ts33'] <- 'LakeSouthNuptse'
HD_met_df$evsites <- evsites

HD_met_df_long <- melt(HD_met_df)

metorder_index <- ggplot(HD_met_df_long, aes(x=sampleID, y=value, fill=evsites)) +
  facet_wrap(~variable, ncol=1, scales='free') +
  geom_point(size=4, shape=21) +
  scale_fill_brewer(palette='Paired') +
  xlab('Sample ID') +
  ylab('Index value') +
  labs(fill='Sites') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        panel.grid.major=element_line(size=0.5, colour='seashell'),
        panel.grid.minor=element_line(size=0.5, colour='seashell'),
        legend.position='bottom',
        legend.text=element_text(size=10)) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

ggsave(plot=metorder_index, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/shansimp_metabar.jpg',
       height=8, width=8, dpi=500)

metotu_df <- data.frame('OTUrankcounts'=rowSums(t(datawide2)))
metotu_df$taxID <- row.names(metotu_df)
dim(metotu_df) #20

metorder_otu <- ggplot(metotu_df, aes(x=taxID, y=OTUrankcounts)) +
  geom_bar(stat='identity') +
  xlab('Taxon ID') +
  ylab('Order OTU counts per taxon ID') +
  ggtitle('Order metabarcoding OTUs') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10))
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabar_shansimp_sampleXorder.jpg',
     height=8, width=15, units='in', res=500)
plot_grid(metorder_index, metorder_otu, ncol=2, align='v', labels=c('a)', 'b)'))
dev.off() 

# ----- MAIN FIGS: by site individually: USE THIS -----
newdf <- as.data.frame.matrix(t(datawide2))
newdf$Order <- row.names(newdf)
newdf$index_for_fig <- c(1:nrow(newdf))
row.names(newdf) <- NULL
metOrderdf <- melt(newdf, id.vars=c('Order', 'index_for_fig'))
metcount_thresh <- 1

write.csv(newdf[,c(20,21,1:19)], 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabar_OrderOTUrelabund_key_data.csv', row.names=FALSE)
length(unique(metOrderdf[metOrderdf$value >=1,]$index_for_fig))# there are 20 Orders with >=1 counts (after changed to >75% ID instead of 70%)
metfig_orders <- unique(metOrderdf[metOrderdf$value >=1,c(1:2)])
metfig_orders[order(metfig_orders$index_for_fig),]

metaboveEBC <- ggplot(metOrderdf[(metOrderdf$variable=='ts13' | metOrderdf$variable=='ts14') & metOrderdf$value >= metcount_thresh, ], 
                   aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Above Everest Base Camp \nTS13 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts13' & metOrderdf$value >= metcount_thresh, ]),
                ' TS14 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts14' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                             legend.position='right')
metKPL1 <- ggplot(metOrderdf[(metOrderdf$variable=='ts15' | metOrderdf$variable=='ts17') & metOrderdf$value >= metcount_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kala Pattar Lake 1 \nTS15 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts15' & metOrderdf$value >= metcount_thresh, ]),
                ' TS17 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts17' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                        legend.position='right')
metKPL2 <- ggplot(metOrderdf[(metOrderdf$variable=='ts18' | metOrderdf$variable=='ts19') & metOrderdf$value >= metcount_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kala Pattar Lake 2 \nTS18 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts18' & metOrderdf$value >= metcount_thresh, ]),
                ' TS19 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts19' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                        legend.position='right')
metKLL1in <- ggplot(metOrderdf[(metOrderdf$variable=='ts20' | metOrderdf$variable=='ts21') & metOrderdf$value >= metcount_thresh, ], 
                 aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 1 (inlet) \nTS20 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts20' & metOrderdf$value >= metcount_thresh, ]),
                ' TS21 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts21' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                                legend.position='right')
metKLL1out <- ggplot(metOrderdf[(metOrderdf$variable=='ts22' | metOrderdf$variable=='ts23') & metOrderdf$value >= metcount_thresh, ], 
                  aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 1 (outlet) \nTS22 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts22' & metOrderdf$value >= metcount_thresh, ]),
                ' TS23 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts23' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                                 legend.position='right')
metKLL2 <- ggplot(metOrderdf[(metOrderdf$variable=='ts24' | metOrderdf$variable=='ts25') & metOrderdf$value >= metcount_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 2 \nTS24 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts24' & metOrderdf$value >= metcount_thresh, ]),
                ' TS25 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts25' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                      legend.position='right')
metKLL3 <- ggplot(metOrderdf[(metOrderdf$variable=='ts26' | metOrderdf$variable=='ts27') & metOrderdf$value >= metcount_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake 3 \nTS26 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts26' & metOrderdf$value >= metcount_thresh, ]),
                ' TS27 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts27' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                      legend.position='right')
metNGMS <- ggplot(metOrderdf[(metOrderdf$variable=='ts28' | metOrderdf$variable=='ts29') & metOrderdf$value >= metcount_thresh, ], 
               aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Nuptse Glacier Mtn Stream \nTS28 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts28' & metOrderdf$value >= metcount_thresh, ]),
                ' TS29 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts29' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                               legend.position='right')
metKLLMS <- ggplot(metOrderdf[(metOrderdf$variable=='ts30' | metOrderdf$variable=='ts31') & metOrderdf$value >= metcount_thresh, ], 
                aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.5) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Kongma La Lake Mtn Stream \nTS30 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts30' & metOrderdf$value >= metcount_thresh, ]),
                ' TS31 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts31' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
                                               legend.position='right')
metLSN <- ggplot(metOrderdf[metOrderdf$variable=='ts33' & metOrderdf$value >= metcount_thresh, ], 
                 aes(x=as.character(index_for_fig), y=value, fill=variable)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width=0.2) +
  xlab(' ') + ylab(' ') + labs(fill='Sample ID') +
  ggtitle(paste('Lake South Nuptse \nTS32 Orders n=', 
                nrow(metOrderdf[metOrderdf$variable=='ts32' & metOrderdf$value >= metcount_thresh, ]),
                ' TS33 Orders n=',
                nrow(metOrderdf[metOrderdf$variable=='ts33' & metOrderdf$value >= metcount_thresh, ]),
                sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15),
        legend.position='right')

jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metOrderOTUall_bysite1-5.jpg',
     height=15, width=15, units='in', res=500)
plot_grid(metaboveEBC, metKPL1, metKPL2, metKLL1in, metKLL1out, ncol=1, align='h')
dev.off()
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metOrderOTUall_bysite6-10.jpg',
     height=15, width=15, units='in', res=500)
plot_grid(metKLL2, metKLL3, metNGMS, metKLLMS, metLSN, ncol=1, align='h')
dev.off()

#try a different layout
panela <- plot_grid(metaboveEBC, metKPL1, metKPL2, ncol=1)
panelb <- plot_grid(metKLL1in, metKLL1out, metKLL2, ncol=1)
panelc <- plot_grid(metKLL3, metNGMS, metKLLMS, ncol=1)
paneld <- plot_grid(metorder_index, metLSN, ncol=1, rel_heights = c(2,1))
jpeg('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metbarOrder_shansimp_OTUall_bysite.jpg',
     height=15, width=25, units='in', res=500)
plot_grid(panela, panelb, panelc, paneld, ncol=4)
dev.off()

# saving individual plots
ggsave(plot=metaboveEBC, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_aboveEBC.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKPL1, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KPL1.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKPL2, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KPL2.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKLL1in, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KLL1in.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKLL1out, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KLL1out.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKLL2, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KLL2.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKLL3, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KLL3.jpg', height=6, width=10, dpi=500)
ggsave(plot=metNGMS, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_NGMS.jpg', height=6, width=10, dpi=500)
ggsave(plot=metKLLMS, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_KLLMS.jpg', height=6, width=10, dpi=500)
ggsave(plot=metLSN, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/metabarOrders_LSN.jpg', height=6, width=10, dpi=500)

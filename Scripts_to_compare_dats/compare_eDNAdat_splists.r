# compare taxa (regardless of coverage, just presence/absence)
# between metabarcoding and WGS dats

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(reshape2)

# ---- read in the datasets -----
metabar_sp <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/metabar_blastresults_onlymitocontigs.csv', header=TRUE)
wgs_sp <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_proxyreadcov_counts.csv', header=TRUE)
spades_sp <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_spadescontigcovs.csv', header=TRUE)

# ---- format data ----
metabar_sp$datatype <- 'metabar_contigs'
wgs_sp$datatype <- 'wgs_reads'
spades_sp$datatype <- 'wgs_contigs'

# set datatype colors, so they are same for all plots
metabar_sp$ptcol <- '#a1d76a'
wgs_sp$ptcol <- '#c51b7d'
spades_sp$ptcol <- '#fee08b'

names(metabar_sp)
names(wgs_sp)
names(spades_sp)

# give same name to rbind(), but remember these are DIFFERENT calcs
colnames(metabar_sp)[4] <- 'conf'
colnames(wgs_sp)[2] <- 'conf'
colnames(spades_sp)[23] <- 'conf'

# give same name to rbind(), but remember these are AVERAGES
colnames(wgs_sp)[4] <- 'per_id'
colnames(wgs_sp)[6] <- 'aln_len'

# combine only fullspname, sampleID, dataset type, phylum, class, order columns
# adding back queryID - mostly for the metabarcoding data, so it's possible to track
# back to the contig and contig sequence info, which are not kept in these taxon tables
# the wgs_sp does not keep queryID due to how its confidence score is calculated (averaged over reads with same sp ID 
# 'coverage' is really not correct to measure for those results, even the proxy confidence score is not
# ideal. since we're presenting the results as present vs. not detected instead,
# i'm going to create a dummy column for the wgs_sp dataframe so that i can still rbind the dataframes)
wgs_sp$queryID <- 'dummycolumn'
alldatasets <- do.call(rbind, list(metabar_sp[,c(1:3,4,6,7,19:23)], wgs_sp[,c(1:2,13,3:4,6,8:12)], spades_sp[,c(1,3,5,7,8,20:25)]))

# add locality/site info to aggregate data by site
evsites <- rep('na', nrow(alldatasets))
evsites[alldatasets$sampleID=='ts13' |
          alldatasets$sampleID=='ts14'] <- 'AboveEBC'
evsites[alldatasets$sampleID=='ts15' |
          alldatasets$sampleID=='ts16' |
          alldatasets$sampleID=='ts17'] <- 'KalaPattarLake1'
evsites[alldatasets$sampleID=='ts18' |
          alldatasets$sampleID=='ts19'] <- 'KalaPattarLake2'
evsites[alldatasets$sampleID=='ts20' |
          alldatasets$sampleID=='ts21'] <- 'KongmaLaLake1_inlet'
evsites[alldatasets$sampleID=='ts22' |
          alldatasets$sampleID=='ts23'] <- 'KongmaLaLake1_outlet'
evsites[alldatasets$sampleID=='ts24' |
          alldatasets$sampleID=='ts25'] <- 'KongmaLaLake2'
evsites[alldatasets$sampleID=='ts26' |
          alldatasets$sampleID=='ts27'] <- 'KongmaLaLake3'
evsites[alldatasets$sampleID=='ts28' |
          alldatasets$sampleID=='ts29'] <- 'NuptseGlacierMtnStream'
evsites[alldatasets$sampleID=='ts30' |
          alldatasets$sampleID=='ts31'] <- 'KongmaLaLakeMtnStream'
evsites[alldatasets$sampleID=='ts32' |
          alldatasets$sampleID=='ts33'] <- 'LakeSouthNuptse'
evsites[alldatasets$sampleID=='ts35'] <- 'moss'
alldatasets$evsites <- evsites

# now, same dot plots as with the separate datasets
# ts32 keeps being plotted out of order, going to set the order:
alldatasets$sampleIDorder <- factor(alldatasets$sampleID, 
                               levels=c('ts13', 'ts14', 'ts15', 'ts16', 'ts17',
                                        'ts18', 'ts19', 'ts20', 'ts21', 'ts22',
                                        'ts23', 'ts24', 'ts25', 'ts26', 'ts27',
                                        'ts28', 'ts29', 'ts30', 'ts31', 'ts32',
                                        'ts33', 'ts35'))
mypal <- unique(alldatasets[,c('datatype', 'ptcol')])
palval <- mypal$ptcol
names(palval) = mypal$datatype

# making percent similarity threshold >90% to remove lower matches
per_id_thresh = 75

# ------- SUPP FIGS: by site -----
# taking out TS16 for Kala Pattar Lake 1 site
# taking out TS35
# taking out bacteria categories
# IMPORTANT!!! FILTERING BASED ON TAXA PRESENT WITH >70% ID or whatever set abofe in per_id_thresh variable
# now trying with >75% (all wgs hits are 80%+ so this only removes some of the poor matches for metabarcoding data)
forsites_df <- alldatasets[alldatasets$sampleID != 'ts16' &
                             alldatasets$sampleID != 'ts35' &
                             alldatasets$phylum != 'Proteobacteria' &
                             alldatasets$phylum != 'Bacteroidetes' &
                             alldatasets$phylum != 'uncultured bacterium' &
                             alldatasets$phylum != 'unknown' &
                             alldatasets$phylum != 'unranked' &
                             alldatasets$phylum != 'unclassified amoebozoa' &
                             alldatasets$phylum != 'unclassified Arthropoda' &
                             alldatasets$phylum != 'unclassified invertebrate' &
                             alldatasets$phylum != 'Choanozoa' & #not enough hits to mean much
                             alldatasets$phylum != 'Ochryophyta' & #cannot find info about this phyla, maybe same or misspelling of Ochrophyta?
                             alldatasets$per_id >= 75 &
                             !is.na(alldatasets$phylum), ]
unique(forsites_df$phylum)
# add Kingdom classifications
kingdom <- rep('na', nrow(forsites_df))
kingdom[forsites_df$phylum=='Streptophyta' |
          forsites_df$phylum=='Rhodophyta' |
          forsites_df$phylum=='Chlorophyta' |
          forsites_df$phylum=='Glaucophyta' |
          forsites_df$phylum=='Cryptophyta'] <- 'Plantae'
kingdom[forsites_df$phylum=='Ochrophyta' |
          forsites_df$phylum=='Heterokonta' |
          forsites_df$phylum=='Bacillariophyta' |
          forsites_df$phylum=='Oomycota' |
          forsites_df$phylum=='Heterokontophyta' |
          forsites_df$phylum=='Haptista' |
          forsites_df$phylum=='Myzozoa' |
          forsites_df$phylum=='Ciliophora'] <- 'Chromista'
kingdom[forsites_df$phylum=='Mucoromycota' |
          forsites_df$phylum=='Blastocladiomycota' |
          forsites_df$phylum=='Basidiomycota' |
          forsites_df$phylum=='Ascomycota'] <- 'Fungi'
kingdom[forsites_df$phylum=='Cnidaria' |
          forsites_df$phylum=='Chordata' |
          # forsites_df$phylum=='Annelida' |
          # forsites_df$phylum=='Nemertea' |
          forsites_df$phylum=='Arthropoda' |
          forsites_df$phylum=='Rotifera' |
          forsites_df$phylum=='Tardigrada'] <- 'Animalia'
kingdom[forsites_df$phylum=='Discosea' |
          forsites_df$phylum=='Amoebozoa'] <- 'Protista'
kingdom[forsites_df$phylum=='Euglenozoa'] <- 'Excavata'
forsites_df$kingdom <- kingdom
unique(forsites_df$kingdom)

write.csv(forsites_df, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/allEvdata_taxinfo_bysite.csv', row.names=FALSE)

ggplot(forsites_df[forsites_df$per_id >=per_id_thresh,], 
       aes(y=evsites, x=phylum, 
           fill=per_id, size=aln_len)) +
  facet_grid(datatype~kingdom, scales='free_x') +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() +
  ylab('Site') +
  ggtitle(paste('eDNA datasets: site x taxon Phylum \n>',per_id_thresh,'% ID', sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'), 
        legend.position='right')
ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/SITE_alldats_above', per_id_thresh, 'perid_byPhy.jpg', sep=''), height=10, width=15, dpi=500)

# metabarcoding and WGS separately
ggplot(forsites_df[forsites_df$per_id >=per_id_thresh &
                     forsites_df$datatype == 'metabar_contigs',], 
       aes(y=evsites, x=phylum, 
           fill=per_id, size=aln_len)) +
  facet_grid(~kingdom, scales='free_x') +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() +
  ylab('Site') +
  xlab('Phylum') +
  labs(fill='% seq sim', size='Alignment len') +
  ggtitle(paste('Metabarcoding: site x Phylum and Kingdom \n>',per_id_thresh,'% ID', sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'), 
        legend.position='right')
ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/SITE_METABARdats_above', per_id_thresh, 'perid_byPhy.jpg', sep=''), 
       height=6, width=12, dpi=500)

ggplot(forsites_df[forsites_df$per_id >=per_id_thresh &
                               forsites_df$datatype != 'metabar_contigs',], 
                 aes(y=evsites, x=phylum, 
                     fill=per_id, size=aln_len)) +
  facet_grid(~kingdom, scales='free_x') +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() +
  ylab('Site') +
  xlab('Phylum') +
  labs(fill='% seq sim', size='Alignment len') +
  ggtitle(paste('WGS: site x Phylum and Kingdom \n>',per_id_thresh,'% ID', sep='')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'), 
        legend.position='right')
ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/SITE_WGSdats_above', per_id_thresh, 'perid_byPhy.jpg', sep=''), 
       height=6, width=12, dpi=500)

# ---- try other site data plots ----
#another way to split out the stats per site by phylum
uniq_sites <- unique(forsites_df$evsites)

# facet by data type
for(i in 1:length(uniq_sites)){
  plt1 <- ggplot(forsites_df[forsites_df$per_id >= per_id_thresh &
                               forsites_df$evsites == uniq_sites[i],], 
         aes(y=per_id, x=phylum)) +
    facet_grid(~datatype, scales='free') +
    geom_point(size=3, alpha=0.5) +
    ggtitle(paste('Site: ', uniq_sites[i], sep='')) +
    theme(axis.text.x=element_blank(), 
          panel.grid.major=element_line(size=0.5, colour='seashell'),
          axis.title.x =element_blank())
  plt2 <- ggplot(forsites_df[forsites_df$per_id >= per_id_thresh &
                               forsites_df$evsites == uniq_sites[i],], 
                 aes(y=aln_len, x=phylum)) +
    geom_point(size=3, alpha=0.5) +
    facet_grid(~datatype, scales='free') +
    theme(axis.text.x=element_blank(), 
          panel.grid.major=element_line(size=0.5, colour='seashell'),
          axis.title.x =element_blank())
  plt3 <- ggplot(forsites_df[forsites_df$per_id >= per_id_thresh &
                               forsites_df$evsites == uniq_sites[i],], 
                 aes(y=conf, x=phylum)) +
    geom_point(size=3, alpha=0.5) +
    facet_grid(~datatype, scales='free') +
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major=element_line(size=0.5, colour='seashell'))
  jpeg(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/', uniq_sites[i], '_statsbyPHY_facets.jpg', sep=''),
       height=10, width=12, units='in', res=500)
  print(plot_grid(plt1, plt2, plt3, ncol=1, align = 'v'))
  dev.off()
}

# ---- MAIN FIG: metabarcoding contigs: presence-only data viz per site -----
x1 <- table(forsites_df[forsites_df$per_id >=per_id_thresh &
                         forsites_df$datatype=='metabar_contigs', 7], #phylum column
            forsites_df[forsites_df$per_id >=per_id_thresh &
                          forsites_df$datatype=='metabar_contigs', 12]) #evsites column
x2 <- melt(x1)

kingdom2 <- rep('none', nrow(x2))
kingdom2[x2$Var1=='Streptophyta' |
           x2$Var1=='Rhodophyta' |
           x2$Var1=='Chlorophyta' |
           x2$Var1=='Glaucophyta' |
           x2$Var1=='Cryptophyta'] <- 'Plantae'
kingdom2[x2$Var1=='Ochrophyta' |
          x2$Var1=='Heterokonta' |
          x2$Var1=='Bacillariophyta' |
          x2$Var1=='Oomycota' |
          x2$Var1=='Heterokontophyta' |
          x2$Var1=='Haptista' |
          x2$Var1=='Myzozoa' |
          x2$Var1=='Ciliophora'] <- 'Chromista'
kingdom2[x2$Var1=='Mucoromycota' |
          x2$Var1=='Blastocladiomycota' |
          x2$Var1=='Basidiomycota' |
          x2$Var1=='Ascomycota'] <- 'Fungi'
kingdom2[x2$Var1=='Cnidaria' |
          x2$Var1=='Chordata' |
          # x2$Var1=='Annelida' |
          # x2$Var1=='Nemertea' |
          x2$Var1=='Arthropoda' |
          x2$Var1=='Rotifera' |
          x2$Var1=='Tardigrada'] <- 'Animalia'
kingdom2[x2$Var1=='Discosea' |
          x2$Var1=='Amoebozoa'] <- 'Protista'
kingdom2[x2$Var1=='Euglenozoa'] <- 'Excavata'
x2$kingdom2 <- kingdom2
unique(x2$kingdom2)
x2_noNA <- x2[x2$kingdom2 != 'none',]
unique(x2_noNA$kingdom2)

# add zeros for 'absence' so there will be empty boxes in the tiles
newval <- c()
for(i in 1:nrow(x2_noNA)){
  therow <- x2_noNA[i,]
  if(therow$value >= 1){
    newval[i] <- '1'
  } else{
    newval[i] <- '0'
  }
}
x3 <- cbind(x2_noNA, newval)
colnames(x3) <- c('phylum', 'evsites', 'value', 'kingdom', 'PresAbs')
unique(x3$kingdom)

# order of phylum by kingdom
x3$phybyking <- factor(x3$phylum, 
                    levels=c('Euglenozoa', 'Chlorophyta','Cryptophyta','Glaucophyta','Rhodophyta','Streptophyta', 
                             'Bacillariophyta','Ciliophora','Haptista','Heterokonta','Heterokontophyta',
                             'Myzozoa','Ochrophyta','Oomycota',
                             'Amoebozoa','Discosea',
                             'Ascomycota','Basidiomycota','Blastocladiomycota','Mucoromycota',
                             'Annelida','Arthropoda','Chordata','Cnidaria','Nemertea','Rotifera','Tardigrada'))

ggplot(x3, aes(x=evsites, y=phybyking)) +
  geom_tile(aes(fill=PresAbs), col='grey') +
  scale_fill_manual(labels=c('Not detected', 'Present'), values = c('white', 'black')) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('Metabarcoding data: Phylum by site') +
  ylab(' ') + xlab('Sites') +
  labs(fill='Status')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/metabar_presenceonly_phylumXsite.jpg',
       height=8, width=6, dpi=500)

# ---- MAIN FIG: WGS contigs and reads: presence-only data viz per site -----
# because i could not get contigs AND reads per ref for all samples
# (either contigs did not form
# or reads were taking too long and never finishing blast parser)
# it would be wrong to compare presence vs. not identified for wgs contigs vs. reads
# it makes more sense to compare presence vs. not identified for wgs vs. metabarcoding
# mostly to drive home the point that different datasets pick up different results
# and might be useful for 1) learning more collectively and 2) driving the direction for future 
# taxa/genes to target
y1 <- table(forsites_df[forsites_df$per_id >=per_id_thresh &
                          (forsites_df$datatype=='wgs_contigs' | forsites_df$datatype=='wgs_reads'), 7], #phylum col
            forsites_df[forsites_df$per_id >=per_id_thresh &
                          (forsites_df$datatype=='wgs_contigs' | forsites_df$datatype=='wgs_reads'), 12]) #evsites col
y2 <- melt(y1)

kingdom2 <- rep('none', nrow(y2))
kingdom2[y2$Var1=='Streptophyta' |
           y2$Var1=='Rhodophyta' |
           y2$Var1=='Chlorophyta' |
           y2$Var1=='Glaucophyta' |
           y2$Var1=='Cryptophyta'] <- 'Plantae'
kingdom2[y2$Var1=='Ochrophyta' |
           y2$Var1=='Heterokonta' |
           y2$Var1=='Bacillariophyta' |
           y2$Var1=='Oomycota' |
           y2$Var1=='Heterokontophyta' |
           y2$Var1=='Haptista' |
           y2$Var1=='Myzozoa' |
           y2$Var1=='Ciliophora'] <- 'Chromista'
kingdom2[y2$Var1=='Mucoromycota' |
           y2$Var1=='Blastocladiomycota' |
           y2$Var1=='Basidiomycota' |
           y2$Var1=='Ascomycota'] <- 'Fungi'
kingdom2[y2$Var1=='Cnidaria' |
           y2$Var1=='Chordata' |
           # y2$Var1=='Annelida' |
           # y2$Var1=='Nemertea' |
           y2$Var1=='Arthropoda' |
           y2$Var1=='Rotifera' |
           y2$Var1=='Tardigrada'] <- 'Animalia'
kingdom2[y2$Var1=='Discosea' |
           y2$Var1=='Amoebozoa'] <- 'Protista'
kingdom2[y2$Var1=='Euglenozoa'] <- 'Excavata'
y2$kingdom2 <- kingdom2
unique(y2$kingdom2)
y2_noNA <- y2[y2$kingdom2 != 'none',]
unique(y2_noNA$kingdom2)

# add zeros for 'absence' so there will be empty boxes in the tiles
newval <- c()
for(i in 1:nrow(y2_noNA)){
  therow <- y2_noNA[i,]
  if(therow$value >= 1){
    newval[i] <- '1'
  } else{
    newval[i] <- '0'
  }
}
y3 <- cbind(y2_noNA, newval)
colnames(y3) <- c('phylum', 'evsites', 'value', 'kingdom', 'PresAbs')
unique(y3$kingdom)

# order of phylum by kingdom
y3$phybyking <- factor(y3$phylum, 
                       levels=c('Euglenozoa', 'Chlorophyta','Cryptophyta','Glaucophyta','Rhodophyta','Streptophyta', 
                                'Bacillariophyta','Ciliophora','Haptista','Heterokonta','Heterokontophyta',
                                'Myzozoa','Ochrophyta','Oomycota',
                                'Amoebozoa','Discosea',
                                'Ascomycota','Basidiomycota','Blastocladiomycota','Mucoromycota',
                                'Annelida','Arthropoda','Chordata','Cnidaria','Nemertea','Rotifera','Tardigrada'))

ggplot(y3, aes(x=evsites, y=phybyking)) +
  geom_tile(aes(fill=PresAbs), col='grey') +
  scale_fill_manual(labels=c('Not detected', 'Present'), values = c('white', 'black')) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('WGS data: Phylum by site') +
  ylab(' ') + xlab('Sites') +
  labs(fill='Status')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/wgscontigread_presenceonly_phylumXsite.jpg',
       height=8, width=6, dpi=500)


# ---- metabarcoding contigs: species hits ----
# 'species' => 97% seq sim
# subset contigs that have >97% ID
# tally contig read totals (coverage) total
# and per phylum and kingdom
forreadtally <- forsites_df[forsites_df$per_id >=97 & forsites_df$datatype=='metabar_contigs', ]
View(forreadtally)

# ---- WGS contigs and reads: species hits ----
# 'species' => 97% seq sim
forreadtally_wgs <- forsites_df[forsites_df$per_id >=97 &
                               (forsites_df$datatype=='wgs_contigs' | forsites_df$datatype=='wgs_reads'), ]
View(forreadtally_wgs)
wgs_uniqsites <- unique(forreadtally_wgs$evsites)
# for(i in 1:length(wgs_uniqsites)){
for(i in c(1:5,8,9)){
  thedat1 <- forreadtally_wgs[forreadtally_wgs$kingdom != 'Plantae' & 
                                forreadtally_wgs$evsites == wgs_uniqsites[i],]
  thedat2 <- forreadtally_wgs[forreadtally_wgs$kingdom == 'Plantae' &
                                forreadtally_wgs$evsites == wgs_uniqsites[i],]
  plt1 <- ggplot(thedat1, 
                 aes(x=fullspname, y=per_id, fill=conf)) +
    facet_grid(~kingdom, scales='free') +
    geom_point(shape=21, alpha=0.8, size=3) +
    scale_fill_viridis_c() +
    ggtitle(paste('Site ', wgs_uniqsites[i], '\nWGS species >97%ID', sep='')) +
    xlab(' ') +
    ylab('Sequence similarity') +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
  plt2 <- ggplot(thedat2, 
                 aes(x=fullspname, y=per_id, fill=conf)) +
    facet_grid(~kingdom) +
    geom_point(shape=21, alpha=0.8, size=3) +
    scale_fill_viridis_c() +
    xlab('Species name') +
    ylab('Sequence similarity') +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
  jpeg(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/wgs_', wgs_uniqsites[i],'_sphitXkingdom.jpg', sep=''),
       height=8, width=8, units='in', res=500)
  print(plot_grid(plt1, plt2, ncol=1))
  dev.off()
}
# i=6 and 7 don't have Plantae
# i=10 doesn't have Chromista or Animalia, just Plantae
for(i in c(6,7)){
  thedat1 <- forreadtally_wgs[forreadtally_wgs$kingdom != 'Plantae' & 
                                forreadtally_wgs$evsites == wgs_uniqsites[i],]
  plt1 <- ggplot(thedat1, 
                 aes(x=fullspname, y=per_id, fill=conf)) +
    facet_grid(~kingdom, scales='free') +
    geom_point(shape=21, alpha=0.8, size=3) +
    scale_fill_viridis_c() +
    ggtitle(paste('Site ', wgs_uniqsites[i], '\nWGS species >97%ID', sep='')) +
    xlab('Species name') +
    ylab('Sequence similarity') +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
  jpeg(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/wgs_', wgs_uniqsites[i],'_sphitXkingdom.jpg', sep=''),
       height=8, width=8, units='in', res=500)
  print(plot_grid(plt1, ncol=1))
  dev.off()
}

thedat_10 <- forreadtally_wgs[forreadtally_wgs$kingdom == 'Plantae' &
                              forreadtally_wgs$evsites == 'KongmaLaLake3',]
plt1 <- ggplot(thedat_10, 
               aes(x=fullspname, y=per_id, fill=conf)) +
  facet_grid(~kingdom, scales='free') +
  geom_point(shape=21, alpha=0.8, size=3) +
  scale_fill_viridis_c() +
  ggtitle(paste('Site KongmaLaLake3 \nWGS species >97%ID', sep='')) +
  xlab('Species name') +
  ylab('Sequence similarity') +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
jpeg(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/eDNA_site_statplts/wgs_KongmaLaLake3_sphitXkingdom.jpg', sep=''),
     height=8, width=8, units='in', res=500)
plot_grid(plt1, ncol=1)
dev.off()


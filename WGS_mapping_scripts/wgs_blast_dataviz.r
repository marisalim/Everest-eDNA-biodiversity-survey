## analyze and plot WGS metagenomic eDNA data
## spades contigs (only a few) or reads (merged + unmerged) blast outputs

library(taxize)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)
library(RColorBrewer)

# ---- Load input files ----
# grab all files per ref and combine
multmerge <- function(mypath){
  filenames <- list.files(path=mypath, pattern='_parsed.txt', full.names=TRUE)
  datalist <- lapply(filenames, function(x){read.csv(file=x, header=TRUE, sep='\t', row.names=NULL)})
  Reduce(function(x,y) {rbind(x,y)}, datalist)
}

chiro_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/chiro_map_results/')
epio_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/epio_map_results/')
forfi_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/forfi_map_results/')
gallus_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/gallus_map_results/')
hypsib_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/hypsib_map_results/')
taky_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/taky_map_results/')
tetra_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/tetra_map_results/')
wadi_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wadi_map_results/')
roti_dat <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/roti_map_results/')
juni_dat <-  multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/juni_map_results/')

# juni_dat just want the chloroplast hits
juni_dat2 <- juni_dat[grep('chloroplast', juni_dat$stitle), ]

# then combine all into 1 file
dflist <- list(chiro_dat, epio_dat, forfi_dat, gallus_dat, hypsib_dat, 
               taky_dat, tetra_dat, wadi_dat, roti_dat, juni_dat2)
allref_dat <- do.call(rbind, dflist)
dim(allref_dat)

# ---- filter out any matches to human/primates ----
# what human matches are there? how many?
allref_dat$fullspname <- paste(allref_dat$genus, allref_dat$species, sep=' ')

allref_dat_nohuman <- allref_dat[allref_dat$genus != 'Homo' & allref_dat$species != 'sapiens' &
                                   allref_dat$species != 'Gorilla' & allref_dat$species != 'Hylobates' &
                                   allref_dat$genus != 'Pongo' & allref_dat$species != 'Pan' &
                                   allref_dat$genus != 'Human' & allref_dat$fullspname != 'Pan troglodytes' &
                                   allref_dat$fullspname != 'Callithrix jacchus' & allref_dat$fullspname != 'MACACA MULATTA' &
                                   allref_dat$fullspname != 'Chlorocebus aethiops' & allref_dat$fullspname != 'Macaca mulatta' &
                                   allref_dat$fullspname != 'Rhesus Macaque' & allref_dat$fullspname != 'PREDICTED: Homo' &
                                   allref_dat$fullspname != 'Canis Familiaris,', ]
dim(allref_dat_nohuman)

# ---- taxize ----
uniq_sp <- unique(allref_dat_nohuman$fullspname)
length(uniq_sp) #now 1276
# # only need to run once! unless you're updating taxonomy
thetaxdf <- data.frame()
# for(i in 1:length(uniq_sp)){
#   thetax <- tax_name(query=uniq_sp[i], get=c('phylum', 'class', 'order'), db='ncbi')
#   thetaxdf <- rbind(thetaxdf, thetax)
#   Sys.sleep(5) #rate limiting bc NCBI ENTREZ only allows 10 requests/second
# }
dim(thetaxdf)

# save so don't need to re-run!!
write.csv(thetaxdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_sphits_taxdf.csv', row.names=FALSE)

# TO DO just for roti new hits
# have to fix the taxa NAs - but do not write over FIRSTVER! otherwise have to
# redo all the ones you already fixed! delete duplicate entries
# just run for roti taxa
# roti <- allref_dat_nohuman[allref_dat_nohuman$ref == 'roti', ]
# uniq_sp <- unique(roti$fullspname)
# length(uniq_sp) #now 187
# thetaxdf <- data.frame()
# for(i in 1:length(uniq_sp)){
#   thetax <- tax_name(query=uniq_sp[i], get=c('phylum', 'class', 'order'), db='ncbi')
#   thetaxdf <- rbind(thetaxdf, thetax)
#   Sys.sleep(5) #rate limiting bc NCBI ENTREZ only allows 10 requests/second
# }
# write.csv(thetaxdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/ROTI_mapped_sphits_taxdf.csv', row.names=FALSE)

# TO DO just for juni chloroplast new hits
# have to fix the taxa NAs - but do not write over FIRSTVER! otherwise have to
# redo all the ones you already fixed! delete duplicate entries
# juni <- allref_dat_nohuman[allref_dat_nohuman$ref == 'juni', ]
# uniq_sp <- unique(juni$fullspname)
# length(uniq_sp) #now 222
# thetaxdf <- data.frame()
# for(i in 1:length(uniq_sp)){
#   thetax <- tax_name(query=uniq_sp[i], get=c('phylum', 'class', 'order'), db='ncbi')
#   thetaxdf <- rbind(thetaxdf, thetax)
#   Sys.sleep(5) #rate limiting bc NCBI ENTREZ only allows 10 requests/second
# }
# write.csv(thetaxdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/JUNI_mapped_sphits_taxdf.csv', row.names=FALSE)


# ---- filter taxize results ---- 
# open above csv and fix NAs in taxize names
mapped_taxa <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_sphits_taxdf_FINAL.csv', header=TRUE)
colnames(mapped_taxa) <- c('db', 'fullspname', 'phylum', 'class', 'order')

allref_dat_taxondat <- merge(x=allref_dat_nohuman, y=mapped_taxa, by.x='fullspname')

# save file!
write.csv(allref_dat_taxondat, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_FULLblasttaxtabs.csv', row.names=FALSE)

# --- *START HERE* ----
# run from here if already ran above code!
### run this read.csv() line if you have not already run above code to generate the table!
allref_dat_taxondat <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_FULLblasttaxtabs.csv', header=TRUE)
above200 <- allref_dat_taxondat[allref_dat_taxondat$aln_len >= 200,]

# ---- Tally sp hits by taxon group ----
allref_dat_taxondat_phylum <- as.data.frame.matrix(table(allref_dat_taxondat$phylum, allref_dat_taxondat$sampleID))
allref_dat_taxondat_phylumdf <- data.frame('Phylum'=row.names(allref_dat_taxondat_phylum), allref_dat_taxondat_phylum)
rownames(allref_dat_taxondat_phylumdf) <- NULL
View(allref_dat_taxondat_phylumdf)

allref_dat_taxondat_class <- as.data.frame.matrix(table(allref_dat_taxondat$class, allref_dat_taxondat$sampleID))
allref_dat_taxondat_classdf <- data.frame('class'=row.names(allref_dat_taxondat_class), allref_dat_taxondat_class)
rownames(allref_dat_taxondat_classdf) <- NULL
View(allref_dat_taxondat_classdf)

write.csv(allref_dat_taxondat_phylumdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_PHYLUMXsamptable.csv', row.names=FALSE)
write.csv(allref_dat_taxondat_classdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_CLASSXsamptable.csv', row.names=FALSE)

ggplot(allref_dat_taxondat, aes(x=phylum)) + 
  geom_histogram(stat='count') + 
  theme(axis.text.x=element_text(angle=45, hjust=1))

# ---- >200bp: plot histograms of basic blast outputs ----
#there should be a cut off for confidence threshold...>200bp
ggplot(allref_dat_taxondat, aes(x=sampleID, y=aln_len)) + 
  geom_point(alpha=0.6) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('Range in blast aln len for \nmapped reads or contigs')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_alnlenrange.jpg', height=6, width=6, dpi=500)

ggplot(above200, aes(x=sampleID, y=aln_len)) + 
  geom_point(alpha=0.6) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('Range in blast aln len for \nmapped reads or contigs')

# phylum histogram
ggplot(above200, aes(x=phylum)) + 
  geom_histogram(stat='count') + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('Histogram of phyla counts \n(includes duplicate counts at this point)')

# %id histogram
ggplot(above200, aes(x=sampleID, y=per_id, 
                     shape=querytype, fill=ref)) + 
  geom_point(alpha=0.6) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle('Range in blast perc id for \nmapped reads or contigs')
# ----- ts28 baetis check -----
baetis <- allref_dat_taxondat[allref_dat_taxondat$sampleID == 'ts28' & allref_dat_taxondat$genus == 'Baetis' & allref_dat_taxondat$aln_len >= 200,]
distinct_baetis <- baetis %>%
  distinct(queryID, .keep_all=TRUE)
write.csv(distinct_baetis, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/ts28_distinctBaetisreadsandcontigs.csv', row.names=FALSE)
# check the sequences in geneious
# with de novo assembly - get 4 contigs, primarily 2 contigs of ~1000 and ~390 reads each
# the reads themselves have relatively few differences, but the consensus sequences
# for these contigs are awful and no good species ID
# the highest % id for seq reads is ~96% to Baetis pseudorhodani, which is endemic to places in Spain
# so, these seqs are highly likely Baetis genus, but that's all we can say

# ----- amoeba check -----
amoeba <- allref_dat_taxondat[allref_dat_taxondat$aln_len >= 200 &
                                allref_dat_taxondat$phylum == 'Amoebozoa',]
dim(amoeba)
# 82-86% id, just 1 read per sample for ts15, 20, 32...so low, wouldn't expect to see these in krona plts from Genewiz


# ---- try 1: try to calculate proxy for coverage...meh ----
# per ref, per sample, are there >1 reads/contigs per species hit?
# note: there will likely be overlap between refs 
uniq_ref <- sort(unique(above200$ref))
uniq_samp <- sort(unique(above200$sampleID))

for(i in 1:length(uniq_ref)){
# for(i in 1:2){ ## for testing
  for(j in 1:length(uniq_samp)){
  # for(j in 1:2){ ## for testing
    subdat <- above200[above200$ref == uniq_ref[i] & above200$sampleID == uniq_samp[j], ]
    if(nrow(subdat) >= 1){
      subplt <- ggplot(subdat, aes(x=fullspname)) + 
        geom_dotplot(stackgroups=TRUE, binwidth = 0.5, alpha=0.5, binpositions='all', fill='orange') +
        facet_grid(querytype~phylum, scales='free') + scale_y_continuous(breaks=NULL) +
        ylab('Count per query type') + xlab('Blast sp hit name') +
        ggtitle(paste('BWA ref:', uniq_ref[i], '\nSample:', uniq_samp[j], sep=' ')) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
      ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/map_proxycov_plts/', 
                   uniq_samp[j], '_', uniq_ref[i], '_proxycov.jpg', sep=''),
             height=5, width=6, dpi=500)
    }
  }
}

# ---- try 2: now, try tally consensus across all refs...meh ----

### instead of looping through uniq refs, loop thru uniq sample and check
# consensus across refs if there are duplicate matches
# but same species across refs, does not mean same read/contig - should see multiple dots per ref per sample in this case
uniq_samp <- sort(unique(above200$sampleID))
for(i in 1:length(uniq_samp)){
  subdat <- above200[above200$sampleID == uniq_samp[i], ]
  if(nrow(subdat) >= 1){
    subplt <- ggplot(subdat, aes(x=fullspname, fill=ref)) + 
      geom_dotplot(stackgroups=TRUE, dotsize=0.5,
                   binpositions='all', binwidth=0.5) +
      facet_grid(querytype~phylum, scales='free') + scale_y_continuous(breaks=NULL) +
      scale_fill_brewer(palette='Paired') +
      ylab('Count per query type') + xlab('Blast sp hit name') +
      ggtitle(paste('Sample:', uniq_samp[i], sep=' ')) +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
    ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/map_consensuscheck_plts/', 
                 uniq_samp[i], '_consensuscheck.jpg', sep=''),
           height=8, width=10, dpi=500)
  }
}

# same as above, but only for arthropods (most common)
# now facet by order!
for(i in 1:length(uniq_samp)){
  subdat <- above200[above200$sampleID == uniq_samp[i] & above200$phylum == 'Arthropoda', ]
  if(nrow(subdat) >= 1){
    subplt <- ggplot(subdat, aes(x=fullspname, fill=ref)) + 
      geom_dotplot(stackgroups=TRUE, dotsize=0.4,
                   binpositions='all', binwidth=0.5) +
      facet_grid(querytype~order, scales='free') + scale_y_continuous(breaks=NULL) +
      scale_fill_brewer(palette='Paired') +
      ylab('Count per query type') + xlab('Blast sp hit name') +
      ggtitle(paste('Arthropoda for Sample:', uniq_samp[i], sep=' ')) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/map_consensuscheck_plts/', 
                 uniq_samp[i], '_arthro_consensuscheck.jpg', sep=''),
           height=8, width=10, dpi=500)
  }
}

# ---- try 3: hm, but want to see blast stats ----
# changing ggplot section of loop

uniq_samp <- sort(unique(above200$sampleID))
for(i in 1:length(uniq_samp)){
  subdat <- above200[above200$sampleID == uniq_samp[i], ]
  if(nrow(subdat) >= 1){
    subplt <- ggplot(subdat, aes(x=per_id, y=aln_len, size=qcovs, 
                                   col=ref, shape=querytype)) +
      geom_point() +
      scale_color_brewer(palette='Paired') +
      ylab('Blast alignment length (bp)') + xlab('Blast % seq match') +
      ggtitle(paste('Blast stats for sample:', uniq_samp[i], sep=' ')) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
    ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/map_blaststat_plts/', 
                 uniq_samp[i], '_blaststats.jpg', sep=''),
           height=8, width=10, dpi=500)
  }
}

# ---- try 4a: ok, improving. here's my calculation for proxy read 'coverage' ----
# to remove duplicate counts per read, keep 1st line for unique reads 
# (looks like the read match irrespective of ref used is the same)

# let's try the blast stats only for distinct matches this time
onlydistinctmatches <- above200 %>%
  distinct(genbankID, per_id, aln_len, num_mismatch, evalue, bitscore, qlen, qcovs, .keep_all = TRUE)

uniq_samp <- sort(unique(onlydistinctmatches$sampleID))
finaldat_df <- data.frame()
for(i in 1:length(uniq_samp)){
  subdat <- droplevels(onlydistinctmatches[onlydistinctmatches$sampleID == uniq_samp[i] & onlydistinctmatches$querytype=='read', ])
  if(nrow(subdat) >= 1){
    # per sample, want the sum of same hits now that duplicates removed
    # still need to separate between contigs and reads (otherwise likely double counting)
    # but since contigs have their own cov from spades, just tally read proxy cov here
    proxycovs_df <- data.frame()
    proxycovs <- as.data.frame(table(subdat$fullspname))
    proxycovs$sampleID = uniq_samp[i]
    colnames(proxycovs) <- c('fullspname', 'proxycov', 'sampleID')
    
    # add average and sd of percent ID, aln len for these matches
    uniq_fullsp <- unique(subdat$fullspname)
    avg_sd_df <- data.frame()
    for(k in 1:length(uniq_fullsp)){
      newdat <- subdat[subdat$fullspname == uniq_fullsp[k], ]
      if(nrow(newdat) != 1){
        perid_avg <- mean(newdat$per_id)
        perid_sd <- round(sd(newdat$per_id),3)
        alnlen_avg <- mean(newdat$aln_len)
        alnlen_sd <- round(sd(newdat$aln_len),3)
        newdat2 <- data.frame('fullspname'=uniq_fullsp[k],
                              'Avg_per_id'=perid_avg, 'SD_per_id'=perid_sd, 'Avg_aln_len'=alnlen_avg, 'SD_aln_len'=alnlen_sd,
                              'queryID'=newdat$queryID)
        avg_sd_df <- rbind(avg_sd_df, newdat2)
      }else if(nrow(newdat) == 1){
        newdat3 <- data.frame('fullspname'=uniq_fullsp[k],
                              'Avg_per_id'=newdat$per_id, 'SD_per_id'=999, 'Avg_aln_len'=newdat$aln_len, 'SD_aln_len'=999)
        avg_sd_df <- rbind(avg_sd_df, newdat3)
      }
    }
    # save dfs and output at end of loop..to save as csv
    proxycovs_df <- rbind(proxycovs_df, proxycovs)
    finaldat <- merge(proxycovs_df, avg_sd_df, by='fullspname')
    finaldat_df <- rbind(finaldat_df, finaldat)
    
    # visually, this is what the results look like...
    subdat2 <- droplevels(onlydistinctmatches[onlydistinctmatches$sampleID == uniq_samp[i], ])
    ggplot(subdat2, aes(x=per_id, y=aln_len, col=qcovs,
                                 shape=querytype)) +
      geom_point(alpha=0.8, size=5) +
      scale_color_viridis_c() +
      geom_text(aes(label = fullspname), size = 3,
                color='slategray') +
      ylab('Blast alignment length (bp)') + xlab('Blast % seq match') +
      ggtitle(paste('Blast stats for sample:', uniq_samp[i], sep=' ')) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(paste('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/map_blaststat_plts/',
                   uniq_samp[i], '_blaststats.jpg', sep=''),
             height=8, width=10, dpi=500)
  }
}

head(finaldat_df)
hist(finaldat_df$proxycov) # there's quite an outlier in 'coverage'
finaldat_df[finaldat_df$proxycov >5, ]

#want to add back taxon info...
mapped_taxa <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_sphits_taxdf_FINAL.csv', header=TRUE)
colnames(mapped_taxa) <- c('db', 'fullspname', 'phylum', 'class', 'order')
proxycov_df2 <- merge(x=finaldat_df, y=mapped_taxa[,c(2:5)], by.x='fullspname')
head(proxycov_df2)
# change 999 to NA
proxycov_df2[proxycov_df2=='999'] <- 'NA'

write.csv(proxycov_df2, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_proxyreadcov_counts.csv', row.names=FALSE)

# let's look at phylum x sampleID for WGS sp list (based on read matches)
wgs_sp <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_proxyreadcov_counts.csv', header=TRUE)
ggplot(wgs_sp, aes(y=sampleID, x=phylum, fill=Avg_per_id, size=Avg_aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS read matches by taxon Phylum') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_reads_byPHY.jpg', height=6, width=8, dpi=500)

ggplot(wgs_sp, aes(y=sampleID, x=class, fill=Avg_per_id, size=Avg_aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS read matches by taxon Class') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_reads_byCLASS.jpg', height=6, width=15, dpi=500)

ggplot(wgs_sp, aes(y=sampleID, x=order, fill=Avg_per_id, size=Avg_aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS read matches by taxon Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_reads_byORDER.jpg', height=6, width=20, dpi=500)

# ---- try 4b: now, the contigs have coverage stats based on spades assembly ----
# it's a kmer coverage, and needs to be extracted from the queryID label
# need good way to explain how it is interpreted
# so, now, let's look at the contig labels and isolate the kmer cov per hit

spadescov_df <- data.frame()
uniq_samp <- sort(unique(onlydistinctmatches$sampleID))
for(i in 1:length(uniq_samp)){
  subdat <- onlydistinctmatches[onlydistinctmatches$sampleID == uniq_samp[i] & onlydistinctmatches$querytype == 'contig', ]
  if(nrow(subdat) >= 1){
    subdat$spadescov <- as.data.frame(matrix(unlist(strsplit(x=as.character(subdat$queryID), split='_')), ncol=6, byrow=TRUE))$V6
    spadescov_df <- rbind(spadescov_df, subdat)
  }
}

head(spadescov_df[,c(1,3,5,7:8,23)])
write.csv(spadescov_df, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_spadescontigcovs.csv', row.names=FALSE)

# same plots now for spades contigs
spades_sp <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/mapped_spadescontigcovs.csv', header=TRUE)
ggplot(spades_sp, aes(y=sampleID, x=phylum, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS contigs matches by taxon Phylum') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_contigs_byPHY.jpg', height=6, width=8, dpi=500)

ggplot(spades_sp, aes(y=sampleID, x=class, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS contigs matches by taxon Class') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_contigs_byCLASS.jpg', height=6, width=15, dpi=500)

ggplot(spades_sp, aes(y=sampleID, x=order, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('WGS contigs matches by taxon Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_contigs_byORDER.jpg', height=6, width=20, dpi=500)

# ---- SUPP FIGURE: btw, what taxa did each ref mito genome catch? ----
# i used a variety of refs, mostly focused on inverts
# but became obvious early on that the mappings were grabbing fairly off-target/non-specific sequences as well
# not a bad thing, since i got more per ref! however, also suggests a ref per taxon not necessarily needed
# so, let's look at the taxon spread per ref
# gonna subset this out by the phylum/class/order tags 

# this time, i just want to see taxon breadth per unique ref, not sample
ref_distincttaxa <- above200 %>%
  distinct(fullspname, ref, phylum, class, order, .keep_all = TRUE)

# there are sometimes multiple sp per phylum/class/order
# in this case, the point will have many overlapping and thus is darker
# light points have fewer or only 1 sp hit
# clean up the phylum categories - don't want the unknown/Na/bacteria
ref_distincttaxa_clean <- ref_distincttaxa[ref_distincttaxa$sampleID != 'ts16' &
                                             ref_distincttaxa$sampleID != 'ts35' &
                                             ref_distincttaxa$phylum != 'Proteobacteria' &
                                             ref_distincttaxa$phylum != 'Bacteroidetes' &
                                             ref_distincttaxa$phylum != 'uncultured bacterium' &
                                             ref_distincttaxa$phylum != 'unknown' &
                                             ref_distincttaxa$phylum != 'unranked' &
                                             ref_distincttaxa$phylum != 'unclassified amoebozoa' &
                                             ref_distincttaxa$phylum != 'unclassified Arthropoda' &
                                             ref_distincttaxa$phylum != 'unclassified invertebrate' &
                                             ref_distincttaxa$phylum != 'Choanozoa' & #not enough hits to mean much
                                             ref_distincttaxa$phylum != 'Ochryophyta' & #cannot find info about this phyla, maybe same or misspelling of Ochrophyta?
                             !is.na(ref_distincttaxa$phylum), ]
unique(ref_distincttaxa_clean$phylum)

# edit ref names - full instead of abbrev
ref_distincttaxa_clean$reforder <- factor(ref_distincttaxa_clean$ref, 
                                    levels=c('juni', 'roti', 'hypsib',
                                             'forfi', 'epio', 'chiro',
                                             'wadi', 'taky', 'tetra',
                                             'gallus'))

ref_distincttaxa_clean$phybyking <- factor(ref_distincttaxa_clean$phylum, 
                       levels=c('Euglenozoa', 'Chlorophyta','Cryptophyta','Glaucophyta','Rhodophyta','Streptophyta', 
                                'Bacillariophyta','Ciliophora','Haptista','Heterokontophyta',
                                'Myzozoa','Ochrophyta',
                                'Amoebozoa','Discosea',
                                'Ascomycota','Basidiomycota','Mucoromycota',
                                'Arthropoda','Chordata','Cnidaria','Rotifera','Tardigrada'))

ggplot(ref_distincttaxa_clean, aes(x=reforder, y=phybyking)) +
  # geom_point(size=3, alpha=0.3) +
  geom_point(size=3) +
  ggtitle('all WGS hits per ref by Phylum \n (mitochondrial or chloroplast ref)') +
  xlab('Mapping reference species') +
  ylab(' ') +
  scale_x_discrete(labels=c('Juniperus recurva','Rotaria rotatoria',
                            'Hypsibius dujardini', 'Lithobius forficatus',
                            'Epiophlebia superstes',
                            'Chironomus tepperi', 
                            'Wadicosa fidelis',
                            'Takydromus amurensis', 
                            'Tetraogallus tibetanus',
                            'Gallus gallus')) +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/BWA_mitogenome_mapping/wgs_reftaxonhits_byPHY.jpg', height=6, width=8, dpi=500)

# wadicosa and chironomus picked up most taxonomic breadth
# tetraogallus, takydromus, and hypsibius were more limited
# picked up amoebozoa, fungus, inverts, a few chordates, the odd cnidaria, algaes, rotifer, possible tardigrade









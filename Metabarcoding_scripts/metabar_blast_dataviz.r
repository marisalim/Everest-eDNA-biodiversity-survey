## analyze and plot metabarcoding eDNA data
## geneious contig and kraken read (merged + unmerged) blast outputs

library(taxize)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)

# ---- Load input files ----
# create single df per parameter set
multmerge <- function(mypath){
  filenames <- list.files(path=mypath, pattern='_parsed.txt', full.names=TRUE)
  datalist <- lapply(filenames, function(x){read.csv(file=x, header=TRUE, sep='\t', row.names=NULL)})
  Reduce(function(x,y) {rbind(x,y)}, datalist)
  }

geneiouscontigs_qcov80perid70df <- multmerge('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/qcov80_perid70/')
dim(geneiouscontigs_qcov80perid70df)

# ---- Filter 1: count rows with Homo sapiens matches per sample and remove rows ----
print(paste('Number of contigs with human species IDs:', nrow(geneiouscontigs_qcov80perid70df[geneiouscontigs_qcov80perid70df$genus == 'Homo' & geneiouscontigs_qcov80perid70df$species == 'sapiens' | geneiouscontigs_qcov80perid70df$species == 'Gorilla' | geneiouscontigs_qcov80perid70df$species == 'Hylobates' ,]), sep=' '))
# 20
geneiouscontigs_qcov80perid70df_nohuman <- geneiouscontigs_qcov80perid70df[geneiouscontigs_qcov80perid70df$genus != 'Homo' & geneiouscontigs_qcov80perid70df$species != 'sapiens' & geneiouscontigs_qcov80perid70df$species != 'Gorilla' & geneiouscontigs_qcov80perid70df$species != 'Hylobates', ]
dim(geneiouscontigs_qcov80perid70df_nohuman)

# ---- Filter 2: classify higher order taxonomy ----
# for this, we will use the taxize package to add higher taxon names
# had to set up API key in r environment: ENTREZ_KEY = "89fa37adc3d9f98c653cd7fbc1160b3e7f09"

# since this takes a while, will only run taxize for unique species names and then merge with blast output after
geneiouscontigs_qcov80perid70df_nohuman$fullspname <- paste(geneiouscontigs_qcov80perid70df_nohuman$genus, geneiouscontigs_qcov80perid70df_nohuman$species, sep=' ')
eneiouscontigs_qcov80perid70df_sp <- unique(geneiouscontigs_qcov80perid70df_nohuman$fullspname)
unique_sp <- union(geneiouscontigs_qcov80perid70df_sp, geneiouscontigs_qcov90perid70df_sp)
length(unique_sp)

thetaxdf <- data.frame()
for(i in 1:length(unique_sp)){
  thetax <- tax_name(query=unique_sp[i], get=c('phylum', 'class', 'order'), db='ncbi')
  thetaxdf <- rbind(thetaxdf, thetax)
  Sys.sleep(5) #rate limiting bc NCBI ENTREZ only allows 10 requests/second
}
dim(thetaxdf)

# save so don't need to re-run!!
write.csv(thetaxdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_contig_sphits_taxdf.csv', row.names=FALSE)

# ---- Merge tax ids with the original blast output file ----
# match columns $fullspname and $query (need to make these the same col name to merge)
geneious_taxa <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_contig_sphits_taxdf.csv', header=TRUE)

colnames(geneious_taxa) <- c('db', 'fullspname', 'phylum', 'class', 'order')

geneious_qcov80perid70_taxondat <- merge(x=geneiouscontigs_qcov80perid70df_nohuman, y=geneious_taxa, by.x='fullspname')

# ---- Search for NA taxonomies and fix them! ----
# ------ 1. use _NAs to find all the cases ------
# these rows have NAs
geneious_qcov80perid70_NAphylum <- geneious_qcov80perid70_taxondat[is.na(geneious_qcov80perid70_taxondat$phylum),] 
geneious_qcov80perid70_NAclass <- geneious_qcov80perid70_taxondat[is.na(geneious_qcov80perid70_taxondat$class),] 
geneious_qcov80perid70_NAorder <- geneious_qcov80perid70_taxondat[is.na(geneious_qcov80perid70_taxondat$order),] 

sp_NAphylum <- unique(geneious_qcov80perid70_NAphylum$fullspname)
## some of these could be called 'unclassified bacteria', 'unclassified Arthropoda', 'unknown'
# [1] "Arthropoda environmental"          "Chironominae sp."                  "Curvibacter putative"             
# [4] "Ectocarpus crouaniorum"            "Florenciella parvula"              "Frankineae bacterium"             
# [7] "Hincksia mitchelliae"              "Kuckuckia spinosa"                 "Lagenidium humanum"               
# [10] "Mammalian expression"              "N.winogradskyi DNA"                "Neorhizobium galegae,"            
# [13] "Nitzschia cf."                     "Stella sp."                        "Tardibacter chloracetimidivorans,"
# [16] "TPA_exp: Gynuella"                 "Uncultured bacterium"              "uncultured Sphingopyxis"  

sp_NAclass <- unique(geneious_qcov80perid70_NAclass$fullspname)
# [1] "Arthropoda environmental"          "Chironominae sp."                  "Curvibacter putative"             
# [4] "Frankineae bacterium"              "Mammalian expression"              "N.winogradskyi DNA"               
# [7] "Neorhizobium galegae,"             "Nitzschia cf."                     "Proteobacteria bacterium"         
# [10] "Rhodothermaceae bacterium"         "Stella sp."                        "Tardibacter chloracetimidivorans,"
# [13] "TPA_exp: Gynuella"                 "Uncultured bacterium"              "uncultured Sphingopyxis"    

sp_NAorder <- unique(geneious_qcov80perid70_NAorder$fullspname)
# [1] "Actinobacteria bacterium"          "Arthropoda environmental"          "Betaproteobacteria bacterium"     
# [4] "Chironominae sp."                  "Curvibacter putative"              "Frankineae bacterium"             
# [7] "Fujientomon dicestum"              "Mammalian expression"              "Micavibrio aeruginosavorus"       
# [10] "N.winogradskyi DNA"                "Neorhizobium galegae,"             "Nitzschia cf."                    
# [13] "Phreatobacter cathodiphilus"       "Phreatobacter sp."                 "Phreatobacter stygius"            
# [16] "Polymorphum gilvum"                "Proteobacteria bacterium"          "Stella sp."                       
# [19] "Tardibacter chloracetimidivorans," "TPA_exp: Gynuella"                 "Uncultured bacterium"             
# [22] "uncultured Sphingopyxis"       

# ------ 2. edit taxa files ------
# edit geneious taxa list
# made edits directly in csv file
# some are easily categorized as unclassified..., others had typos or could find the missing taxon label via google search

# ------ 3. redo the merge function ------
geneious_qcov80perid70_taxondat <- merge(x=geneiouscontigs_qcov80perid70df_nohuman, y=geneious_taxa, by.x='fullspname')

# check for NAs
# just run code from step 1 (should get empty dfs because no more NAs!)

# ------ 4. save these tables! ------
write.csv(geneious_qcov80perid70_taxondat, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_qcov80perid70_FULLblasttaxtabs.csv', row.names=FALSE)

# --- *START HERE* ----
# run from here if already ran above code!
### run these read.csv() lines if you have not already run above code to generate the tables!
geneious_qcov80perid70_taxondat <- read.csv('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_qcov80perid70_FULLblasttaxtabs.csv', header=TRUE)

# ---- Tally sp hits by taxon group ----
geneious_qcov80perid70_phylum <- as.data.frame.matrix(table(geneious_qcov80perid70_taxondat$phylum, geneious_qcov80perid70_taxondat$sampleID))
geneious_qcov80perid70_phylumdf <- data.frame('Phylum'=row.names(geneious_qcov80perid70_phylum), geneious_qcov80perid70_phylum)
rownames(geneious_qcov80perid70_phylumdf) <- NULL
geneious_qcov80perid70_phylumdf

geneious_qcov80perid70_class <- as.data.frame.matrix(table(geneious_qcov80perid70_taxondat$class, geneious_qcov80perid70_taxondat$sampleID))
geneious_qcov80perid70_classdf <- data.frame('Class'=row.names(geneious_qcov80perid70_class), geneious_qcov80perid70_class)
rownames(geneious_qcov80perid70_classdf) <- NULL
geneious_qcov80perid70_classdf

write.csv(geneious_qcov80perid70_phylumdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_qcov80perid70_PHYLUMXsamptable.csv', row.names=FALSE)
write.csv(geneious_qcov80perid70_classdf, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/geneious_qcov80perid70_CLASSxsamptable.csv', row.names=FALSE)

# ---- bacteria filter, Tally sp hits by taxon group, count read cov for geneious contigs ----
gene8070_sums <- geneious_qcov80perid70_taxondat %>% 
  group_by(sampleID, phylum) %>% 
  summarise(total = sum(readspercontig))

ggplot(gene8070_sums, aes(x=sampleID, y=total, group=phylum, fill=phylum)) + 
  geom_bar(stat='identity') +
  coord_flip() +
  # geom_text(aes(label = comma(total)), size = 2, position=position_stack(vjust=0.5)) +
  ggtitle('Reads/contig phylum count \n(geneious qcov80perid70)') +
  ylab('Reads per contig') +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=15),
        axis.text.y= element_text(size=15),
        axis.title=element_text(size=20, face='bold'), legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/readspercontig_byPhylumALL_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

#### WITHOUT bacteria
#### keeping mitochondrial hits only
dim(geneious_qcov80perid70_taxondat)
nobacdat <- geneious_qcov80perid70_taxondat[grep('mitochondrial', geneious_qcov80perid70_taxondat$stitle), ]
nobacdat2 <- geneious_qcov80perid70_taxondat[grep('mitochondrion', geneious_qcov80perid70_taxondat$stitle), ]
nobacdat_final <- do.call(rbind, list(nobacdat, nobacdat2))
dim(nobacdat_final)
write.csv(nobacdat_final, 'Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/metabar_blastresults_onlymitocontigs.csv', row.names=FALSE)

gene8070_sums2 <- nobacdat_final %>% 
  group_by(sampleID, phylum, class, order) %>% 
  summarise(total = sum(readspercontig))
ggplot(gene8070_sums2, aes(x=sampleID, y=total, group=phylum, fill=phylum)) + 
  # facet_grid(~phylum, scales='free') +
  geom_bar(stat='identity') +
  coord_flip() +
  # geom_text(aes(label = phylum), size = 2, position=position_stack(vjust=0.5), alpha=0.4) +
  ggtitle('Reads/contig phylum count \n(geneious qcov80perid70)') +
  ylab('Reads per contig') +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=15),
        axis.text.y= element_text(size=15),
        axis.title=element_text(size=20, face='bold'), legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/readspercontig_byPhylumnobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

ggplot(gene8070_sums2, aes(x=sampleID, y=total, group=class, fill=class)) + 
  geom_bar(stat='identity') +
  coord_flip() +
  # geom_text(aes(label = comma(total)), size = 2, position=position_stack(vjust=0.5)) +
  ggtitle('Reads/contig class count \n(geneious qcov80perid70)') +
  ylab('Reads per contig') +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=15),
        axis.text.y= element_text(size=15),
        axis.title=element_text(size=20, face='bold'), legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/readspercontig_byClassnobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

ggplot(gene8070_sums2, aes(x=sampleID, y=total, group=order, fill=order)) + 
  geom_bar(stat='identity') +
  coord_flip() +
  # geom_text(aes(label = comma(total)), size = 2, position=position_stack(vjust=0.5)) +
  ggtitle('Reads/contig order count \n(geneious qcov80perid70)') +
  ylab('Reads per contig') +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=15),
        axis.text.y= element_text(size=15),
        axis.title=element_text(size=20, face='bold'), legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/readspercontig_byOrdernobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

# same plots as WGS reads for phylum/class/order, but now metabarcoding for contigs
ggplot(nobacdat_final, aes(y=sampleID, x=phylum, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('Metabarcoding contig matches by taxon Phylum') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/metabar_contigs_byPHY.jpg', height=6, width=8, dpi=500)

ggplot(nobacdat_final, aes(y=sampleID, x=class, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('Metabarcoding contig matches by taxon Class') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/metabar_contigs_byCLASS.jpg', height=6, width=8, dpi=500)

ggplot(nobacdat_final, aes(y=sampleID, x=order, fill=per_id, size=aln_len)) +
  geom_point(alpha=0.4, pch=21) +
  scale_fill_viridis_c() + ggtitle('Metabarcoding contig matches by taxon Order') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'))
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/metabar_contigs_byORDER.jpg', height=6, width=8, dpi=500)

## want to plot out the coverage too -- but size param scales make it difficult to see
ggplot(nobacdat_final, aes(x=readspercontig, y=sampleID)) +
  geom_point() 
## most are <5000x

# ------ General data plots ------
# histogram of aln_len
gene8070_plt <- ggplot(geneious_qcov80perid70_taxondat) + 
  geom_histogram(aes(x=aln_len), bins=10) + ggtitle('Blast result alignment length \ngeneious qcov80 perID70')

# contigs or reads per sample
plot(table(geneious_qcov80perid70_taxondat$sampleID))

# look at perID, aln_len, qcov per sample
# without Proteobacteria
ggplot(geneious_qcov80perid70_taxondat[geneious_qcov80perid70_taxondat$phylum != 'Proteobacteria', ], aes(x=sampleID, y=per_id/100, col=phylum)) +
  geom_point(alpha=0.7, size=5) +
  # geom_hline(yintercept=c(0.8, 0.9), lty=2, alpha=0.6, col='red') +
  scale_y_continuous(limits=c(0.7,1),
                     labels=scales::percent) +
  ylab('Percent sequence similarity') +
  ggtitle('Range in contig % id by sample, excluding Proteobacteria') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'),
        legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/perID_byPhylumnoprobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

ggplot(geneious_qcov80perid70_taxondat[geneious_qcov80perid70_taxondat$phylum != 'Proteobacteria', ], aes(x=sampleID, y=qcovs/100, col=phylum)) +
  geom_point(alpha=0.7, size=5) +
  scale_y_continuous(limits=c(0.8,1),
                     labels=scales::percent) +
  ylab('Percent query alignment coverage') +
  ggtitle('Range in contig % qcov by sample, excluding Proteobacteria') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'),
        legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/qcov_byPhylumnoprobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

ggplot(geneious_qcov80perid70_taxondat[geneious_qcov80perid70_taxondat$phylum != 'Proteobacteria', ], aes(x=sampleID, y=aln_len, col=phylum)) +
  geom_point(alpha=0.7, size=5) +
  ylab('Blast match alignment length (bp)') +
  ggtitle('Range in contig aln length by sample, excluding Proteobacteria') +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        panel.grid.major=element_line(size=0.5, colour='seashell'),
        legend.position='bottom')
ggsave('Desktop/Everest_metagen_dat/MetagenomicsAnalyses/METABARCODING_analyses/Tracie_metabar_geneious_contigs/alnlen_byPhylumnoprobac_geneiousqcov80perid70.jpg', 
       height=8, width=12, dpi=500)

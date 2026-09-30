library(dada2)
library(ShortRead)
library(phyloseq)

## Quality-filter the sequencing reads and create Amplicon Sequence Variant (ASV) tables with DADA2

path <- "~/Documents/Panama_Hypoxia_TentExperiment/extraction_blank_analysis/cutadapt"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Perform filtering and trimming
filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
# Learn the Error Rates, it TAKES TIME!
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# Inspecting the dada-class object returned by dada:
dadaFs[[1]]
# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
head(track)
#281 tabled reads in extraction blank (Kit-3)

write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)
saveRDS(seqtab.nochim, "~/Documents/Panama_Hypoxia_TentExperiment/extraction_blank_analysis/seqtab.rds") 


## Assign taxonomy in DADA2
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
# FIX the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "silva_otu_table.txt",sep="\t",col.names=NA)

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Tent_Expm_Metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

# remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") 
get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom") 
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")
get_taxa_unique(ps, "Family") 
get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom") 
ps 

# export tables
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"silva_nochloronomito_otu_table-withBLANK.txt",sep="\t",col.names=NA)
write.table(taxon,"silva_nochloronomito_taxa_table-withBLANK.txt",sep="\t",col.names=NA)



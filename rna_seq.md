# RNA-Seq analysis

## Load salmon

```
module load Salmon
```

## Downloading the data

For this tutorial we will use the test data from [this](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004393) paper:

> Malachi Griffith*, Jason R. Walker, Nicholas C. Spies, Benjamin J. Ainscough, Obi L. Griffith*. 2015. Informatics for RNA-seq: A web resource for analysis on the cloud. PLoS Comp Biol. 11(8):e1004393.

The test data consists of two commercially available RNA samples: Universal Human Reference (UHR) and Human Brain Reference (HBR). The UHR is total RNA isolated from a diverse set of 10 cancer cell lines. The HBR is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old.

In addition, a spike-in control was used. Specifically we added an aliquot of the ERCC ExFold RNA Spike-In Control Mixes to each sample. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies). This range allows us to test the degree to which the RNA-seq assay (including all laboratory and analysis steps) accurately reflects the relative abundance of transcript species within a sample. There are two 'mixes' of these transcripts to allow an assessment of differential expression output between samples if you put one mix in each of your two comparisons. In our case, Mix1 was added to the UHR sample, and Mix2 was added to the HBR sample. We also have 3 complete experimental replicates for each sample. This allows us to assess the technical variability of our overall process of producing RNA-seq data in the lab.

For all libraries we prepared low-throughput (Set A) TruSeq Stranded Total RNA Sample Prep Kit libraries with Ribo-Zero Gold to remove both cytoplasmic and mitochondrial rRNA. Triplicate, indexed libraries were made starting with 100ng Agilent/Strategene Universal Human Reference total RNA and 100ng Ambion Human Brain Reference total RNA. The Universal Human Reference replicates received 2 ul of 1:1000 ERCC Mix 1. The Human Brain Reference replicates received 1:1000 ERCC Mix 2. The libraries were quantified with KAPA Library Quantification qPCR and adjusted to the appropriate concentration for sequencing. The triplicate, indexed libraries were then pooled prior to sequencing. Each pool of three replicate libraries were sequenced across 2 lanes of a HiSeq 2000 using paired-end sequence chemistry with 100bp read lengths.

So to summarize we have:

* UHR + ERCC Spike-In Mix1, Replicate 1
* UHR + ERCC Spike-In Mix1, Replicate 2
* UHR + ERCC Spike-In Mix1, Replicate 3
* HBR + ERCC Spike-In Mix2, Replicate 1
* HBR + ERCC Spike-In Mix2, Replicate 2
* HBR + ERCC Spike-In Mix2, Replicate 3

You can download the data from [here](http://139.162.178.46/files/tutorials/toy_rna.tar.gz).

Unpack the data and go into the toy_rna directory

```
tar xzf toy_rna.tar.gz
cd toy_rna
```

## Indexing transcriptome

```
salmon index -t chr22_transcripts.fa -i chr22_index
```

## Quantify reads using salmon

```bash
for i in *_R1.fastq.gz
do
   prefix=$(basename $i _R1.fastq.gz)
   salmon quant -i chr22_index --libType A \
          -1 ${prefix}_R1.fastq.gz -2 ${prefix}_R2.fastq.gz -o quant/${prefix};
done
```

This loop simply goes through each sample and invokes salmon using fairly basic options:

* The -i argument tells salmon where to find the index
* --libType A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
* The -1 and -2 arguments tell salmon where to find the left and right reads for this sample (notice, salmon will accept gzipped FASTQ files directly).
* the -o argument specifies the directory where salmon’s quantification results sould be written.

Salmon exposes many different options to the user that enable extra features or modify default behavior. However, the purpose and behavior of all of those options is beyond the scope of this introductory tutorial. You can read about salmon’s many options in the [documentation](http://salmon.readthedocs.io/en/latest/).

After the salmon commands finish running, you should have a directory named `quant`, which will have a sub-directory for each sample. These sub-directories contain the quantification results of salmon, as well as a lot of other information salmon records about the sample and the run. The main output file (called quant.sf) is rather self-explanatory. For example, take a peek at the quantification file for sample `HBR_Rep1` in `quant/HBR_Rep1/quant.sf` and you’ll see a simple TSV format file listing the name (Name) of each transcript, its length (Length), effective length (EffectiveLength) (more details on this in the documentation), and its abundance in terms of Transcripts Per Million (TPM) and estimated number of reads (NumReads) originating from this transcript.

## Import read counts using tximport

Using the tximport R package, you can import salmon’s transcript-level quantifications and optionally aggregate them to the gene level for gene-level differential expression analysis. 

First, open up your favourite R IDE and install the necessary packages:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("GenomicFeatures")

install.packages("readr")
```

Then load the modules:

```R
library(tximport)
library(GenomicFeatures)
library(readr)
```

Salmon did the quantifiation of the transcript level. We want to see which genes are differentially expressed, so we need to link the transcript names to the gene names. We can use our .gtf annotation for that, and the GenomicFeatures package:

```R
txdb <- makeTxDbFromGFF("chr22_genes.gtf")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)
```

Now we can import the salmon quantification. First, download the file with sample descriptions from [here](https://raw.githubusercontent.com/HadrienG/tutorials/master/data/samples.txt) and put it in the toy_rna directory. Then, use that file to load the corresponding quantification data.

```R
samples <- read.table("samples.txt", header = TRUE)
files <- file.path("quant", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
```

Take a look at the data:

```R
head(txi.salmon$counts)
```

## Differential expression using DESeq2

Install the necessary package:

```R
biocLite('DESeq2')
```

Then load it:

```R
library(DESeq2)
```

Instantiate the DESeqDataSet and generate result table. See `?DESeqDataSetFromTximport` and `?DESeq` for more information about the steps performed by the program.

```R
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
dds <- DESeq(dds)
res <- results(dds)
```

Run the `summary` command to get an idea of how many genes are up- and downregulated between the two conditions:

`summary(res)`

DESeq uses a negative binomial distribution. Such distributions have two parameters: mean and dispersion. The dispersion is a parameter describing how much the variance deviates from the mean.

You can read more about the methods used by DESeq2 in the [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) or the [vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq/inst/doc/DESeq.pdf)

Plot dispersions:

```R
plotDispEsts(dds, main="Dispersion plot")
```

For clustering and heatmaps, we need to log transform our data:

```R
rld <- rlogTransformation(dds)
head(assay(rld))
```

Then, we create a sample distance heatmap:

```R
library(RColorBrewer)
library(gplots) # you may need to install this package

(mycols <- brewer.pal(8, "Dark2")[1:length(unique(samples$condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[samples$condition],
          RowSideColors=mycols[samples$condition],
          margin=c(10, 10), main="Sample Distance Matrix")
```

We can also plot a PCA:

```R
DESeq2::plotPCA(rld, intgroup="condition")
```

It is time to look at some p-values:

```R
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
```

Examine plot of p-values, the MA plot and the Volcano Plot:

```R
hist(res$pvalue, breaks=50, col="grey")
DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)

# Volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
```

## KEGG pathway analysis

As always, install and load the necessary packages:

```R
biocLite("AnnotationDbi")
biocLite("org.Hs.eg.db")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")

library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)
```

Let’s use the `mapIds` function to add more columns to the results. The row.names of our results table has the Ensembl gene ID (our key), so we need to specify  `keytype=ENSEMBL`. The column argument tells the `mapIds` function which information we want, and the `multiVals` argument tells the function what to do if there are multiple possible values for a single input value. Here we ask to just give us back the first one that occurs in the database. Let’s get the Entrez IDs, gene symbols, and full gene names.

```R
res$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

head(res)
```

We’re going to use the [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package for pathway analysis, and the [pathview](https://bioconductor.org/packages/release/bioc/html/pathview.html) package to draw a pathway diagram.

The gageData package has pre-compiled databases mapping genes to KEGG pathways and GO terms for common organisms:

```R
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
```

Run the pathway analysis. See help on the gage function with `?gage`. Specifically, you might want to try changing the value of same.dir.

```R
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)
```

Pull out the top 5 upregulated pathways, then further process that just to get the IDs. We’ll use these KEGG pathway IDs downstream for plotting. The `dplyr` package is required to use the pipe (`%>%`) construct. 

```R
library(dplyr)

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

Finally, the `pathview()` function in the pathview package makes the plots. Let’s write a function so we can loop through and draw plots for the top 5 pathways we created above.

```R
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Unload dplyr since it conflicts with the next line
detach("package:dplyr", unload=T)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
```

#### Thanks

This material was inspired by Stephen Turner's blog post:

> Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview: http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html

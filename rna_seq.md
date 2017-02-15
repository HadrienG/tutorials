# RNA-Seq analysis

## Load salmon

```
module load salmon
```

## Downloading the data.

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

You can download the data from [here](http://139.162.178.46/files/tutorials/toy_rna.tar.gz)

Unpack the data and go into the toy_rna directory

```
tar xzf toy_rna.tar.gz
cd toy_rna
```

## indexing transcriptome

```
salmon index -t chr22_transcripts.fa -i chr22_index
```

## quantify reads using salmon

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

## import read counts using tximport

Using the tximport R package, you can import salmon’s transcript-level quantifications and optionally aggregate them to the gene level for gene-level differential expression analysis

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

Salmon did the quantifiation of the transcript level. We want to see which genes are differentially expressed, so we need to link the transcripts name to the gene names. We can use our .gtf annotation for that, and the GenomicFeatures package:

```R
txdb <- makeTxDbFromGFF("chr22_genes.gtf")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)
```

now we can import the salmon quantification:

```R
samples <- read.table("samples.txt", header = TRUE)

files <- file.path("salmon", samples$quant, "quant.sf")
names(files) <- paste0("sample", 1:6)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
```

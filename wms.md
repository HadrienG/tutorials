# Whole Metagenome Sequencing

In this tutorial you'll analyze a sample from Pig Gut Metagenome.

### Table of Contents

* [Introduction](#introduction)
    * [The Pig Microbiome](#the-pig-microbiome)
    * [Whole Metagenome Sequencing](#whole-metagenome-sequencing)
* [Softwares Required for this Tutorial](#softwares-required-for-this-tutorial)
* [Getting the Data and Checking their Quality](#getting-the-data-and-checking-their-quality)
* [Taxonomic Classification](#taxonomic-classification)
* [Visualization](#visualization)

### Introduction

#### Microbiome used

In this tutorial we will compare samples from the Pig Gut Microbiome to samples from the Human Gut Microbiome. Below you'll find a brief description of the two projects:

> Pig is a main species for livestock and biomedicine. The pig genome sequence was recently reported. To boost research, we established a catalogue of the genes of the gut microbiome based on faecal samples of 287 pigs from France, Denmark and China. More than 7.6 million non-redundant genes representing 719 metagenomic species were identified by deep metagenome sequencing, highlighting more similarities with the human than with the mouse catalogue. The pig and human catalogues share only 12.6 and 9.3 % of their genes, respectively, but 70 and 95% of their functional pathways. The pig gut microbiota is influenced by gender, age and breed. Analysis of the prevalence of antibiotics resistance genes (ARGs) reflected antibiotics supplementation in each farm system, and revealed that non-antibiotics-fed animals still harbour ARGs. The pig catalogue creates a resource for whole metagenomics-based studies, highly valuable for research in biomedicine and for sustainable knowledge-based pig farming

> We are facing a global metabolic health crisis provoked by an obesity epidemic. Here we report the human gut microbial composition in a population sample of 123 non-obese and 169 obese Danish individuals. We find two groups of individuals that differ by the number of gut microbial genes and thus gut bacterial richness. They harbour known and previously unknown bacterial species at different proportions; individuals with a low bacterial richness (23% of the population) are characterized by more marked overall adiposity, insulin resistance and dyslipidaemia and a more pronounced inflammatory phenotype when compared with high bacterial richness individuals. The obese individuals among the former also gain more weight over time. Only a few bacterial species are sufficient to distinguish between individuals with high and low bacterial richness, and even between lean and obese. Our classifications based on variation in the gut microbiome identify subsets of individuals in the general white adult population who may be at increased risk of progressing to adiposity-associated co-morbidities

#### Whole Metagenome Sequencing

Whole Metagenome sequencing (WMS), or shotgun metagenome sequencing, is a relatively new and powerful sequencing approach that provides insight into community biodiversity and function. On the contrary of Metabarcoding, where only a specific region of the bacterial community (the 16s rRNA) is sequenced, WMS aims at sequencing all the genomic material present in the environment.

The choice of shotgun or 16S approaches is usually dictated by the nature of the studies being conducted. For instance, 16S is well suited for analysis of large number of samples, i.e., multiple patients, longitudinal studies, etc. but offers limited taxonomical and functional resolution. WMS is generally more expensive but offers increased resolution, and allows the discovery of archaea and viruses.

### Softwares Required for this Tutorial

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Kraken](https://ccb.jhu.edu/software/kraken/)
<!-- * [scythe](https://github.com/vsbuffalo/scythe) -->
<!-- * kraken_to_krona (script part of [MetLab](https://github.com/norling/metlab))
* [KronaTools](https://github.com/marbl/Krona/wiki) -->

### Getting the Data and Checking their Quality

If you are reading this tutorial online and haven't cloned the directory, first download and unpack the data:

```
wget http://77.235.253.14/metlab/wms.tar
tar xvf wms.tar
cd wms
```

We'll use FastQC to check the quality of our data. FastQC can be downloaded and
ran on a Windows or LINUX computer without installation. It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Start FastQC and select the fastq file you just downloaded with `file -> open`  
What do you think about the quality of the reads? Do they need trimming? Is there still adapters
present? Overrepresented sequences?

Alternatively, run fastqc on the command-line:

`fastqc *.fastq.gz`

If the quality appears to be that good, it's because it was probably the cleaned reads that were deposited into SRA.
We can directly move to the classification step.

### Taxonomic Classification

[Kraken](https://ccb.jhu.edu/software/kraken/) is a system for assigning taxonomic labels to short DNA sequences (i.e. reads)  
Kraken aims to achieve high sensitivity and high speed by utilizing exact alignments of k-mers and a novel classification algorithm (sic).

In short, kraken uses a new approach with exact k-mer matching to assign taxonomy to short reads. It is *extremely* fast compared to traditional
approaches (i.e. Blast)

By default, the authors of kraken built their database based on RefSeq Bacteria, Archea and Viruses. We'll use it for the purpose of this tutorial.

**NOTE: The database may have been installed already! Ask your instructor!**

```bash
# You might not need this step (example if you're working on Uppmax!)
wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
tar xzf minikraken.tgz
$KRAKEN_DB=minikraken_20141208
```

Now run kraken on the reads

```bash
mkdir kraken_results
for i in *_1.fastq.gz
do
    prefix=$(basename $i _1.fastq.gz)
    kraken --db $KRAKEN_DB --threads 2 --fastq-input --gzip-compressed \
        ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz > kraken_results/${prefix}.tab
    kraken-report --db $KRAKEN_DB \
        kraken_results/${prefix}.tab > kraken_results/${prefix}_tax.txt
done
```

which produces a tab-delimited file with an assigned TaxID for each read.

Kraken includes a script called `kraken-report` to transform this file into a "tree" view with the percentage of reads assigned to each taxa. We've run this script at each step in the loop. Take a look at the `_tax.txt` files!

### Abundance estimation using Bracken

Bracken (Bayesian Reestimation of Abundance with KrakEN) is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample

Before starting, you need to install Bracken:

```bash
cd
git clone https://github.com/jenniferlu717/Bracken.git
chmod 755 Bracken/*.py
chmod 755 Bracken/*.pl
export PATH=$PATH:$HOME/Bracken
```

Unfortunately, Uppmax lacks some perl packages necessary for Bracken to work:

Follow the tutorial [here](http://www.uppmax.uu.se/support/faq/software-faq/installing-local-perl-packages/) to install `cpanm`

then install the two perl libraries that are missing:

```bash
cpanm Parallel::ForkManager
cpanm List::MoreUtils
```

Three steps are necessary to set up Kraken abundance estimation.

1. Classify all reads using Kraken and Generate a Kraken report file. We've done this!

2. Search all library input sequences against the database and compute the classifications for each perfect read of ${READ_LENGTH} base pairs from one of the input sequences.


```bash
find -L $KRAKEN_DB/library -name "*.fna" -o -name "*.fa" -o -name "*.fasta" > genomes.list
cat $(grep -v '^#' genomes.list) > genomes.fasta
kraken --db=${KRAKEN_DB} --fasta-input --threads=10 kraken.fasta > database.kraken
count-kmer-abundances.pl --db=${KRAKEN_DB} --read-length=100 database.kraken > database100mers.kraken_cnts
```

3. Generate the kmer distribution file

```bash
python generate_kmer_distribution.py -i database100mers.kraken_cnts -o KMER_DISTR.TXT
```

Now, given the expected kmer distribution for genomes in a kraken database along
with a kraken report file, the number of reads belonging to each species (or
genus) is estimated using the estimate_abundance.py file, run with the
following command line:

`python estimate_abundance.py -i KRAKEN.REPORT -k KMER_DISTR.TXT -o OUTPUT_FILE.TXT`

Run this command for the six `_tax.txt` files that you generated with kraken!

The following required parameters must be specified:
- KRAKEN.REPORT     :: the kraken report generated for a given dataset
- KMER_DISTR.TXT    :: the file generated by generate_kmer_distribution.py
- OUTPUT_FILE.TXT   :: the desired name of the output file to be generated by the code

### Visualization

#### Alternative 1: Pavian

Pavian is a web application for exploring metagenomics classification results.

Install and run Pavian:

(In R or Rstudio)

```R
## Installs required packages from CRAN and Bioconductor
source("https://raw.githubusercontent.com/fbreitwieser/pavian/master/inst/shinyapp/install-pavian.R")
pavian::runApp(port=5000)
```

Pavian will be available at http://127.0.0.1:5000 .

# Whole Metagenome Sequencin

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

The Pig Microbiome:

> Pig is a main species for livestock and biomedicine. The pig genome sequence was recently reported. To boost research, we established a catalogue of the genes of the gut microbiome based on faecal samples of 287 pigs from France, Denmark and China. More than 7.6 million non-redundant genes representing 719 metagenomic species were identified by deep metagenome sequencing, highlighting more similarities with the human than with the mouse catalogue. The pig and human catalogues share only 12.6 and 9.3 % of their genes, respectively, but 70 and 95% of their functional pathways. The pig gut microbiota is influenced by gender, age and breed. Analysis of the prevalence of antibiotics resistance genes (ARGs) reflected antibiotics supplementation in each farm system, and revealed that non-antibiotics-fed animals still harbour ARGs. The pig catalogue creates a resource for whole metagenomics-based studies, highly valuable for research in biomedicine and for sustainable knowledge-based pig farming

The Human Microbiome:

> We are facing a global metabolic health crisis provoked by an obesity epidemic. Here we report the human gut microbial composition in a population sample of 123 non-obese and 169 obese Danish individuals. We find two groups of individuals that differ by the number of gut microbial genes and thus gut bacterial richness. They harbour known and previously unknown bacterial species at different proportions; individuals with a low bacterial richness (23% of the population) are characterized by more marked overall adiposity, insulin resistance and dyslipidaemia and a more pronounced inflammatory phenotype when compared with high bacterial richness individuals. The obese individuals among the former also gain more weight over time. Only a few bacterial species are sufficient to distinguish between individuals with high and low bacterial richness, and even between lean and obese. Our classifications based on variation in the gut microbiome identify subsets of individuals in the general white adult population who may be at increased risk of progressing to adiposity-associated co-morbidities

#### Whole Metagenome Sequencing

Whole Metagenome sequencing (WMS), or shotgun metagenome sequencing, is a relatively new and powerful sequencing approach that provides insight into community biodiversity and function. On the contrary of Metabarcoding, where only a specific region of the bacterial community (the 16s rRNA) is sequenced, WMS aims at sequencing all the genomic material present in the environment.

The choice of shotgun or 16S approaches is usually dictated by the nature of the studies being conducted. For instance, 16S is well suited for analysis of large number of samples, i.e., multiple patients, longitudinal studies, etc. but offers limited taxonomical and functional resolution. WMS is generally more expensive but offers increased resolution, and allows the discovery of viruses as well as other mobile genetic elements.

### Softwares Required for this Tutorial

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Kraken](https://ccb.jhu.edu/software/kraken/)
* [R](https://www.r-project.org/)
* [Pavian](https://github.com/fbreitwieser/pavian)

### Prepare and organise your working directory

```
mkdir wms
cd wms
mkdir data
mkdir results
mkdir scripts
```

### Getting the Data and Checking their Quality

If you are reading this tutorial online and haven't cloned the directory, first download and unpack the data:

```
cd data
wget http://77.235.253.14/metlab/wms.tar
tar xvf wms.tar
cd wms
```

We'll use FastQC to check the quality of our data. FastQC can be downloaded and
run on a Windows or Linux computer without installation. It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Start FastQC and select the fastq file you just downloaded with `file -> open`  
What do you think about the quality of the reads? Do they need trimming? Are there still adapters
present? Overrepresented sequences?

Alternatively, run fastqc on the command-line:

`fastqc *.fastq.gz`

If the quality appears to be good, it's because it was probably the cleaned reads that were deposited into SRA.
We can directly move to the classification step.

### Taxonomic Classification

[Kraken](https://ccb.jhu.edu/software/kraken/) is a system for assigning taxonomic labels to short DNA sequences (i.e. reads)  
Kraken aims to achieve high sensitivity and high speed by utilizing exact alignments of k-mers and a novel classification algorithm (sic).

In short, kraken uses a new approach with exact k-mer matching to assign taxonomy to short reads. It is *extremely* fast compared to traditional
approaches (i.e. BLAST).

By default, the authors of kraken built their database based on RefSeq Bacteria, Archaea and Viruses. We'll use it for the purpose of this tutorial.
We will download a shrinked database (minikraken) provided by Kraken developers that is only 4GB.

```bash
# First we create a databases directory in our home
cd /mnt
sudo mkdir databases
cd databases
# Then we download the minikraken database
sudo wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz
sudo tar xzf minikraken_20171019_4GB.tgz
KRAKEN_DB=/mnt/databases/minikraken_20171013_4GB
cd
```

Now run kraken on the reads

```bash
# In the data/wms directory
cd wms/data/wms
for i in *_1.fastq.gz
do
    prefix=$(basename $i _1.fastq.gz)
    # print which sample is being processed
    echo $prefix
    kraken --db $KRAKEN_DB --threads 2 --fastq-input --gzip-compressed \
        ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz > /home/student/wms/results/${prefix}.tab
    kraken-report --db $KRAKEN_DB \
        /home/student/wms/results/${prefix}.tab > /home/student/wms/results/${prefix}_tax.txt
done
```

which produces a tab-delimited file with an assigned TaxID for each read.

Kraken includes a script called `kraken-report` to transform this file into a "tree" view with the percentage of reads assigned to each taxa. We've run this script at each step in the loop. Take a look at the `_tax.txt` files!

### Visualization with Pavian

Pavian is a web application for exploring metagenomics classification results.

First, go in Rstudio server by typing the address to your server in your browser:

`http://MY_IP_ADDRESS:8787/`

where you replace `MY_IP_ADDRESS` by the IP address of your Virtual Machine.

!!! note
    To access Rstudio server on the virtual machine, you'll need a password
    Ask your instructor for the password!

!!! note
    If you wish, you may work on Rstudio on your own laptop if it is powerful enough.
    You will need an up-to-date version of R, and can install the necessary packages using [this script](https://osf.io/a7kqz/download)


Install and run Pavian:


```R
options(repos = c(CRAN = "http://cran.rstudio.com"))
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```

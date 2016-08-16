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

#### The Pig Microbiome

> Pig is a main species for livestock and biomedicine. The pig genome sequence was recently reported. To boost research, we established a catalogue of the genes of the gut microbiome based on faecal samples of 287 pigs from France, Denmark and China. More than 7.6 million non-redundant genes representing 719 metagenomic species were identified by deep metagenome sequencing, highlighting more similarities with the human than with the mouse catalogue. The pig and human catalogues share only 12.6 and 9.3 % of their genes, respectively, but 70 and 95% of their functional pathways. The pig gut microbiota is influenced by gender, age and breed. Analysis of the prevalence of antibiotics resistance genes (ARGs) reflected antibiotics supplementation in each farm system, and revealed that non-antibiotics-fed animals still harbour ARGs. The pig catalogue creates a resource for whole metagenomics-based studies, highly valuable for research in biomedicine and for sustainable knowledge-based pig farming

To speed up the analysis, we'll only use the first 60K reads from the first sample of the study. The full samples are accessible under BioProject [PRJEB11755](http://www.ncbi.nlm.nih.gov/bioproject/308698)

#### Whole Metagenome Sequencing

Whole Metagenome sequencing (WMS), or shotgun metagenome sequencing, is a relatively new and powerful sequencing approach that provides insight into community biodiversity and function. On the contrary of Metabarcoding, where only a specific region of the bacterial community (the 16s rRNA) is sequenced, WMS aims at sequencing all the genomic material present in the environment.

The choice of shotgun or 16S approaches is usually dictated by the nature of the studies being conducted. For instance, 16S is well suited for analysis of large number of samples, i.e., multiple patients, longitudinal studies, etc. but offers limited taxonomical and functional resolution. WMS is generally more expensive but offers increased resolution, and allows the discovery of archaea and viruses.

### Softwares Required for this Tutorial

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Kraken](https://ccb.jhu.edu/software/kraken/)
* kraken_to_krona (script part of [MetLab](https://github.com/norling/metlab))
* [KronaTools](https://github.com/marbl/Krona/wiki)

### Getting the Data and Checking their Quality

If you are reading this tutorial online and haven't cloned the directory, first download and unzip the data:

```
wget https://github.com/HadrienG/tutorials/raw/master/data/DHN_Pig_60K.fastq.gz
gzip -d DHN_Pig_60K.fastq.gz
```

We'll use FastQC to check the quality of our data. FastQC can be downloaded and
ran on a Windows or LINUX computer without installation. It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Start FastQC and select the fastq file you just downloaded with `file -> open`  
What do you think about the quality of the reads? Do they need trimming? Is there still adapters
present? Overrepresented sequences?

If the quality appears to be that good, it's because it was probably the cleaned reads that were deposited into SRA.
We can directly move to the classification step.

### Taxonomic Classification

[Kraken](https://ccb.jhu.edu/software/kraken/) is a system for assigning taxonomic labels to short DNA sequences (i.e. reads)  
Kraken aims to achieve high sensitivity and high speed by utilizing exact alignments of k-mers and a novel classification algorithm (sic).

In short, kraken uses a new approach with exact k-mer matching to assign taxonomy to short reads. It is *extremely* fast compared to traditional
approaches (i.e. Blast)

By default, the authors of kraken built their database based on RefSeq Bacteria, Archea and Viruses. We'll use it for the purpose of this tutorial.

```
wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz
tar xzf minikraken.tgz
```

Now run kraken on the reads

`kraken --db minikraken_20141208/ --threads 8 --fastq-input DHN_Pig_60K.fastq > DHN_Pig.tab`

which produces a tab-delimited file with an assigned TaxID for each read. Kraken includes a script called `kraken-report` to transform this file into a "tree" view with the percentage of reads assigned to each taxa.

`kraken-report --db minikraken_20141208/ DHN_Pig.tab > DHN_Pig_tax.txt`

Open this file and take a look!

### Visualization

We'll visualize the composition of our datasets using Krona.

Get the script to transform the kraken results in a format Krona can understand

```
wget https://raw.githubusercontent.com/norling/metlab/master/pipeline_scripts/kraken_to_krona.py
chmod 755 kraken_to_krona.py
```

Run the script and Krona

```
./kraken_to_krona.py DHN_Pig_tax.txt
ktImportText DHN_Pig_tax.krona.in
```

And open the generated html file in your browser! What are the most abundant taxa?

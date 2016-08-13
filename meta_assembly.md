# Metagenome assembly

In this tutorial you'll learn how to inspect the quality of High-throughput sequencing and
perform a metagenomic assembly.

We will use data under the accession SRS018585 in the Sequence Read Archive. this sample is
"a Human Metagenome sample from G_DNA_Anterior nares of a male participant in the dbGaP study
HMP Core Microbiome Sampling Protocol A (HMP-A)"

### Table of contents

* [Softwares required for this tutorial](#softwares-required-for-this-tutorial)
* [getting the data](#getting-the-data)
* [quality control](#quality-control)
* [assembly](#assembly)
* [taxonomic classification and visualization](#taxonomic-classification-and-visualization)

### Softwares required for this tutorial

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [sickle](https://github.com/najoshi/sickle)
* [SPAdes](http://bioinf.spbau.ru/en/spades)
* [Blast](https://blast.ncbi.nlm.nih.gov)
* [blobtools](https://drl.github.io/blobtools/)

### getting the data

```
wget http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2
tar xjf SRS018585.tar.bz2
cd SRS018585
```

### quality control

we'll use FastQC to check the quality of our data. FastQC can be downloaded and
ran on a Windows or LINUX computer without installation. It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Start FastQC and select the fastq files you just downloaded with `file -> open`

What is the average read length? The average quality?

Now we'll trim the reads using sickle

```
sickle pe -f SRS018585.denovo_duplicates_marked.trimmed.1.fastq \
-r SRS018585.denovo_duplicates_marked.trimmed.2.fastq -t sanger \
-o SRS018585_trimmed_1.fastq -p SRS018585_trimmed_2.fastq -s unpaired.fastq
```

sickle normally gives you a summary of how many reads were trimmed.

### assembly

SPAdes will be used for the assembly. Since version 3.7, SPAdes includes a metagenomic version of its algorithm, callable
with the option --meta

```
spades.py --meta -1 SRS018585_trimmed_1.fastq -2 SRS018585_trimmed_2.fastq -t 8 -o assembly
```

the resulting assenmbly can be found under assembly/scaffolds.fasta. How many contigs does this assembly contain?
How long is the longest contig and to what organism does it belong to?

### taxonomic classification and visualization

For the vizualisation of the assembly we will use a tool called blobtools.
Blobtools produces "Taxon annotated GC-coverage plots" (TAGC) and was orignially made for
the visualisation of (draft) genome assemblies.  

```
mkdir blobtools && cd $_
blastn -num_threads 8 -db nt -query ../assembly/scaffolds.fasta -out blastresults.txt -outfmt '6 qseqid staxids bitscore'
```

This blast step is necessary to obtain the taxonomic information of your contigs.
It might take a while. Be patient!

```
blobtools create -i ../assembly/scaffolds.fasta -y spades -t blastresults.txt \
    --nodes /export/databases/taxonomy/nodes.dmp \
    --names /export/databases/taxonomy/names.dmp \
    -o scaffolds --title SRS018585
blobtools plot -i scaffolds.blob.BlobDB.json -o scaffolds --title -r family
```

Inspect the plot, what is the most abundant families? try to play with the parameters
(especially `-r`)

# Metagenome assembly and binning

In this tutorial you'll learn how to inspect assemble metagenomic data and retrieve draft genomes from assembled metagenomes

We'll use a mock community of 20 bacteria sequenced using the Illumina HiSeq.
In reality the data were simulated using [InSilicoSeq](http://insilicoseq.readthedocs.io).

The 20 bacteria in the dataset were selected from the [Tara Ocean study](http://ocean-microbiome.embl.de/companion.html) that recovered 957 distinct Metagenome-assembled-genomes (or MAGs) that were previsouly unknown! (full list on [figshare](https://figshare.com/articles/TARA-NON-REDUNDANT-MAGs/4902923/1) )

## Getting the Data

```bash
mkdir -p ~/data
cd ~/data
curl -O -J -L https://osf.io/th9z6/download
curl -O -J -L https://osf.io/k6vme/download
chmod -w tara_reads_R*
```

## Quality Control

we'll use FastQC to check the quality of our data, as well as sickle for trimming the bad quality part of the reads.
If you need a refresher on how and why to check the quality of sequence data, please check the [Quality Control and Trimming](qc) tutorial

```bash
mkdir -p ~/results
cd ~/results
ln -s ~/data/tara_reads_* .
fastqc tara_reads_*.fastq.gz
```

!!! question
    What is the average read length? The average quality?

!!! question
    Compared to single genome sequencing, what graphs differ?


Now we'll trim the reads using sickle

```
sickle pe -f tara_reads_R1.fastq.gz -r tara_reads_R2.fastq.gz -t sanger \
    -o tara_trimmed_R1.fastq -p tara_trimmed_R2.fastq -s /dev/null
```

!!! question
    How many reads were trimmed?

## Assembly

Megahit will be used for the assembly.

```
megahit -1 tara_trimmed_R1.fastq -2 tara_trimmed_R2.fastq -o tara_assembly
```

the resulting assenmbly can be found under `tara_assembly/final.contigs.fa`.

!!! question
    How many contigs does this assembly contain?

## Binning

First we need to map the reads back against the assembly to get coverage information

```bash
ln -s tara_assembly/final.contigs.fa .
bowtie2-build final.contigs.fa final.contigs
bowtie2 -x final.contigs -1 tara_reads_R1.fastq.gz -2 tara_reads_R2.fastq.gz | \
    samtools view -bS -o tara_to_sort.bam
samtools sort tara_to_sort.bam -o tara.bam
samtools index tara.bam
```

then we run metabat

```bash
runMetaBat.sh -m 1500 final.contigs.fa tara.bam
mv final.contigs.fa.metabat-bins1500 metabat
```

!!! question
    How many bins did we obtain?

## Checking the quality of the bins

The first time you run `checkm` you have to create the database

```bash
sudo checkm data setRoot ~/.local/data/checkm
```

```bash
checkm lineage_wf -x fa metabat checkm/
checkm bin_qa_plot -x fa checkm metabat plots
```

!!! question
    Which bins should we keep for downstream analysis?

!!! note
    checkm can plot a lot of metrics. If you have time, check the manual
    and try to produce different plots


!!! warning
    if checkm fails at the phylogeny step, it is likely that your vm doesn't have enough RAM.
    pplacer requires about 35G of RAM to place the bins in the tree of life.

    In that case, execute the following

    `cd ~/results`  
    `curl -O -J -L https://osf.io/xuzhn/download`  
    `tar xzf checkm.tar.gz`  
    `checkm qa checkm/lineage.ms checkm`  

then plot the completeness

```bash
checkm bin_qa_plot -x fa checkm metabat plots
```

and take a look at `plots/bin_qa_plot.png`

## Further reading

* [Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life](https://www.nature.com/articles/s41564-017-0012-7)
* [The reconstruction of 2,631 draft metagenome-assembled genomes from the global oceans](https://www.nature.com/articles/sdata2017203)

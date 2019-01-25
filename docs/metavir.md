# Viral Metagenome from a dolphin sample: hunting for a disease causing virus

In this tutorial you will learn how to investigate metagenomic data and retrieve draft genome from an assembled metagenome.

We will use a real dataset published in 2017 in a study in dolphins,
where fecal samples where prepared for viral metagenomics study.
The dolphin had a self-limiting gastroenteritis of suspected viral origin.

## Getting the Data

First, create an appropriate directory to put the data:
```bash
mkdir -p ~/dolphin/data
cd ~/dolphin/data
```

You can download them from here:
```
curl -O -J -L https://osf.io/4x6qs/download
curl -O -J -L https://osf.io/z2xed/download
```
Alternatively, your instructor will let you know where to get the dataset from.

You should get 2 compressed files:
```
Dol1_S19_L001_R1_001.fastq.gz
Dol1_S19_L001_R2_001.fastq.gz
```

## Quality Control

We will use FastQC to check the quality of our data, as well as sickle for trimming the bad quality part of the reads.
If you need a refresher on how and why to check the quality of sequence data, please check the [Quality Control and Trimming](qc) tutorial

```bash
mkdir -p ~/dolphin/results
cd ~/dolphin/results
ln -s ~/dolphin/data/Dol1* .
fastqc Dol1_*.fastq.gz
```

!!! question
    What is the average read length? The average quality?

!!! question
    Compared to single genome sequencing, which graphs differ?

Removing adapters with Scythe

The first thing we need is the adapters to trim off

```bash
curl -O -J -L https://osf.io/v24pt/download
```

Now we run scythe on both our read files
```bash
scythe -a adapters.fasta -o Dol1_adapt_R1.fastq Dol1_S19_L001_R1_001.fastq.gz
scythe -a adapters.fasta -o Dol1_adapt_R2.fastq Dol1_S19_L001_R2_001.fastq.gz
```

Now we'll trim the reads using sickle

```
sickle pe -f Dol1_adapt_R1.fastq -r Dol1_adapt_R2.fastq -t sanger \
    -o Dol1_trimmed_R1.fastq -p Dol1_trimmed_R2.fastq -s /dev/null -q 25
```

!!! question
    How many reads were trimmed?

## Removing the host sequences by mapping/aligning on the dolphin genome

For this we will use Bowtie2.
We have downloaded the genome of Tursiops truncatus from Ensembl (fasta file).

Alternatively, you can get it from ensembl FTP:
```
curl -O -J -L ftp://ftp.ensembl.org/pub/release-94/fasta/tursiops_truncatus/dna/Tursiops_truncatus.turTru1.dna.toplevel.fa.gz
```

Then we have built the Bowtie2 indexes for this genomes using the following command (do not run it, we have pre-calculated the results for you):
```
bowtie2-build Tursiops_truncatus.turTru1.dna.toplevel.fa Tursiops_truncatus
```

First we will copy the bowtie indexed of the dolphin genome into our results directory:

```
cp /opt/dolphin/Tursiops_truncatus*.bt* .
```

Now we are ready to map our sequencing reads on the dolphin genome:
```
bowtie2
```


## Taxonomic classification of the trimmed reads

We will use Kaiju for the classification of the produced contigs. As we are mainly interested in detecting the viral sequences in our dataset and we want to reduce the computing time and the memory needed, we have built a viruses-only database.

```bash
kaiju -t /opt/kaijudb/nodes.dmp -f /opt/kaijudb/kaiju_db.fmi -i Dol1_trimmed_R1.fastq -j Dol1_trimmed_R2.fastq -o Dol1_trimmed_reads_kaiju.out
```

In order to visualise the results, we will produce a Krona chart

```bash
kaiju2krona -t /opt/kaijudb/nodes.dmp -n /opt/kaijudb/names.dmp -i Dol1_trimmed_reads_kaiju.out -o Dol1_trimmed_reads_kaiju.krona -u

ktImportText -o Dol1_trimmed_reads_kaiju.krona.html Dol1_trimmed_reads_kaiju.krona

```

!!! note
    If we were interested in bacteria and had a databases containing them, we could also use the command kaijuReport to get a text summary. Unfortunately, it does not provide taxonomic levels for viruses.

Then we copy it locally to visualise in our web browser.
From the terminal on our laptop, we move to the place where we want to have the chart (cd) and then copy:
```bash
scp studentX@ebiokit_ip:~/dolphin/results/*krona.html .
```

!!! note
    You need to change the X by your student number and *ebiokit_ip* by the correct ip address provided by your instructor

## Assembly

Megahit will be used for the *de novo* assembly of the metagenome.

```
megahit -1 Dol1_trimmed_R1.fastq -2 Dol1_trimmed_R2.fastq -o Dol1_assembly
```

The resulting assembly can be found under `Dol1_assembly/final.contigs.fa`.

!!! question
    How many contigs does this assembly contain? Is there any long contig?

## Taxonomic classification

We will use Kaiju again with the same viruses-only database for the classification of the produced contigs.

```bash
cd Dol1_assembly

kaiju -t /opt/kaiju/kaijudb/nodes.dmp -f /opt/kaiju/kaijudb/kaiju_db.fmi -i final.contigs.fa -o Dol1_contigs_kaiju.out
```
Then we produce the Krona chart:
```bash
kaiju2krona -t /opt/kaijudb/nodes.dmp -n /opt/kaijudb/names.dmp -i Dol1_contigs_kaiju.out -o Dol1_contigs_kaiju.krona -u

ktImportText -o Dol1_contigs_kaiju.krona.html Dol1_contigs_kaiju.krona
```

!!! question
    Does the classification of contigs produce different results than the classification of reads?

## Extraction of the contig of interest

Let's go for a little practice of your Unix skills!

!!! question
    Find a way to save all the contigs headers. They look like this:
```
>k141_1 flag=1 multi=1.0000 len=301
>k141_2 flag=1 multi=1.0000 len=303
```

**hint 1**: all lines containing '>'

**hint 2**: use grep and the file redirection

Once we have this file, we want to sort all the sequences headers by the sequence length (len=X):

```
cat final.contigs.fa | grep ">" | sed s/len=// | sort -k4,4n | cut -d ' ' -f1 | cut -c2- | tail -1
```

!!! question
    What is the size of the longest contig?

Now that you have identified the sequence header or id of the longest contig, you want to save it to a fasta file.

```
grep -i '>k141_XXX' -A 1 > longest_contig.fasta
```

!!! note
    You need to replace the XXX by the correct header.
    The option -A 1 of grep allows to print 1 line additionally to the matching line (which enables to print the full sequence that corresponds to one line).
    Test with -A 2 (without the redirection to longest_contig.fasta) to see what happens.


Now that you have identified the longest contig, you will check in Kaiju results what was the taxon assigned to this contig.

Have a look at the file **Dol1_contigs_kaiju.out**. It is structured in 3 columns:
classification status (C/U), sequence id, assigned TaxID.

!!! question
    Identify the TaxID of the longest contig and search on NCBI Taxonomy database to which species it corresponds to.



## Genome annotation of the contig of interest

Once the contig to annotate is extracted and saved in the file longest_contig.fasta, we will use Prodigal to detect ORFs (Open Reading Frames) in order to predict genes and their resulting proteins.

```bash
prodigal -i longest_contig.fasta -g 1 -o longest_contig_prodigal.gff -f gff -a longest_contig_prodigal_prot.faa -s longest_contig_prodigal_genes.fasta
```

!!! question
    How many genes and proteins were predicted?

Now that we have done the structural annotation, *i.e.* prediction of the genes genomic location, we will proceed to the functional annotation, *i.e.* predict the proteins functions. We will use Blast to do this.

```
blastp ...
```

## Visualization and manual curation.

If there is some time left, you can download Ugene from the eBioKit and install is on your computer.
Then you can copy the produced gff file and add the functions to the predicted proteins.

## Alternative solution for annotation of prokaryotic and viral genomes

For doing the structural and functional annotation by using one single tool, we recommend to use [Prokka]().

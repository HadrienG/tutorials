# Viral Metagenome from a dolphin sample: hunting for a disease causing virus

In this tutorial you will learn how to investigate metagenomics data and retrieve draft genome from an assembled metagenome.

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

We will use FastQC to check the quality of our data, as well as fastp for trimming the bad quality part of the reads.
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

## Quality control with Fastp
We will removing the adapters and trim by quality.

Now we run fastp our read files
```bash
fastp -i Dol1_S19_L001_R1_001.fastq.gz -o Dol1_trimmed_R1.fastq \
 -I Dol1_S19_L001_R2_001.fastq.gz -O Dol1_trimmed_R2.fastq \
 --detect_adapter_for_pe --length_required 30 \
 --cut_front --cut_tail --cut_mean_quality 10
```

Check the html report produced.

!!! question
    How many reads were trimmed?

## Removing the host sequences by mapping/aligning on the dolphin genome

For this we will use Bowtie2.
We have downloaded the genome of Tursiops truncatus from Ensembl (fasta file). Then we have run the following command to produce the indexes of the dolphin genome for Bowtie2 (do not run it, we have pre-calculated the results for you):

```
bowtie2-build Tursiops_truncatus.turTru1.dna.toplevel.fa Tursiops_truncatus
```

Because this step takes a while, we have precomputed the index files, you can get them from here:

```
curl -O -J -L https://osf.io/wfk9t/download
```

First we will extract the bowtie indexes of the dolphin genome into our results directory:

```
tar -xzvf host_genome.tar.gz
```

Now we are ready to map our sequencing reads on the dolphin genome:
```
bowtie2 -x host_genome/Tursiops_truncatus \
-1 Dol1_trimmed_R1.fastq -2 Dol1_trimmed_R2.fastq \
-S dol_map.sam --un-conc Dol_reads_unmapped.fastq --threads 4
```
!!! question
    How many reads mapped on the dolphin genome?

## Taxonomic classification of the trimmed reads

We will use Kaiju for the classification of the produced contigs. As we are mainly interested in detecting the viral sequences in our dataset and we want to reduce the computing time and the memory needed, we will build a viruses-only database.

```bash
# Kaiju needs to be installed
cd
mkdir -p databases/kaijudb
cd databases/kaijudb
makeDB.sh -v
```

Once the database built, Kaiju tells us that we need only 3 files:
Then other files and directories can be removed

```bash
# Removing un-needed things
rm -r genomes/
rm kaiju_db.bwt
rm kaiju_db.sa
rm merged.dmp
rm kaiju_db.faa
```

We are now ready to run Kaiju on our trimmed reads

```bash
cd dolphin/results
kaiju -t ~/databases/kaijudb/nodes.dmp -f ~/databases/kaijudb/kaiju_db.fmi \
 -i Dol_reads_unmapped.1.fastq -j Dol_reads_unmapped.2.fastq \
 -o Dol1_reads_kaiju.out
```

In order to visualise the results, we will produce a Krona chart.
This step requires to have KronaTools installed:
```
conda install -c bioconda krona
```

```bash
kaiju2krona -t ~/databases/kaijudb/nodes.dmp -n ~/databases/kaijudb/names.dmp \
-i Dol1_reads_kaiju.out -o Dol1_reads_kaiju.krona -u

ktImportText -o Dol1_reads_kaiju.krona.html Dol1_reads_kaiju.krona

```

!!! note
    If we were interested in bacteria and had a databases containing them, we could also use the command kaijuReport to get a text summary. Unfortunately, it does not provide taxonomic levels for viruses.

Then we copy the produced html file locally to visualise in our web browser.


## Assembly

Megahit will be used for the *de novo* assembly of the metagenome.

```
megahit -1 Dol_reads_unmapped.1.fastq -2 Dol_reads_unmapped.2.fastq -o assembly
```

The resulting assembly can be found under `assembly/final.contigs.fa`.

!!! question
    How many contigs does this assembly contain? Is there any long contig?

## Taxonomic classification of contigs

We will use Kaiju again with the same viruses-only database for the classification of the produced contigs.

```bash
cd assembly

kaiju -t ~/databases/kaijudb/nodes.dmp -f ~/databases/kaijudb/kaiju_db.fmi \
 -i final.contigs.fa -o Dol1_contigs_kaiju.out
```
Then we produce the Krona chart:
```bash
kaiju2krona -t ~/databases/kaijudb/nodes.dmp -n ~/databases/kaijudb/names.dmp \
 -i Dol1_contigs_kaiju.out -o Dol1_contigs_kaiju.krona -u

ktImportText -o Dol1_contigs_kaiju.krona.html Dol1_contigs_kaiju.krona
```

!!! question
    Does the classification of contigs produce different results than the classification of reads?

## Extraction of the contig of interest

Let's go for a little practice of your Unix skills!

!!! question
    Find a way to to find the longest contig. They look like this:
```
>k141_1 flag=1 multi=1.0000 len=301
>k141_2 flag=1 multi=1.0000 len=303
```

**hint 1**: all lines containing '>'

**hint 2**: use grep and sed

Once we have this file, we want to sort all the sequences headers by the sequence length (len=X):

```bash
cat final.contigs.fa | grep ">" | sed s/len=// | sort -k4n | tail -1
```

!!! question
    What is the size of the longest contig?

Now that you have identified the sequence header or id of the longest contig, you want to save it to a fasta file.

```bash
grep -i '>k141_XXX' -A 1 final.contigs.fa > longest_contig.fasta
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

Once the contig to annotate is extracted and saved in the file longest_contig.fasta, we will use Prokka to detect ORFs (Open Reading Frames) in order to predict genes and their resulting proteins.

First, go to Uniprot database and retrieve a set of protein sequences
belonging to adenoviruses. Save the file as `adenovirus.faa` and copy it in
your results directory.

```bash
prokka --outdir annotation --kingdom Viruses \
--proteins adenovirus.faa longest_contig.fasta
```

!!! question
    How many genes and proteins were predicted?

## Visualization and manual curation.

If there is some time left, you can visualise the produced annotation (gff file)
in Ugene or Artemis for example.


## Go further with proteins functions

For the predicted proteins that are left "hypotetical", you can try running
[Interproscan](https://www.ebi.ac.uk/interpro/search/sequence-search) on them to get more information on domains and motifs.

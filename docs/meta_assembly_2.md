# Metagenome assembly and binning (continued)

In the previous tutorial  we have seen how to recover draft genomes from raw metagenomic data.
But we do not know what kind of organisms we have yet!

Identyfying an unknown species from a genome is not always easy, especially if the species has not been described before.
In this tutorial, you'll learn a few ways of identyfying (or attempting to identidy) an unkown genome

## Getting the Data

This tutorial focuses on one particular genome bin but you can download the necessary data by executing the code block below.

```bash
mkdir -p ~/mag/data
cd ~/mag/data
# TODO link here
chmod -w *.fa.gz
```

## Ribosomal RNA

The simplest thing we can do is search our draft genome for rRNA genes.
Since these genes are usually quite conserved across species/genera, it could give us a broad idea of our organism

```bash
barrnap -o bin2_rrna.fa bin.2.fa
cat bin2_rrna.fa
```

Now, Use the online Blast service to search similar sequences to the rRNA we obtained.

## Assigning taxonomy to each contig

We'll use diamond against the swissprot database for quickly assigning taxonomy to our contigs.

First, we download and build the database

```bash
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
diamond makedb --in uniprot_sprot.fasta.gz --db uniprot_sprot -p 4
```

then we run diamond

```bash
diamond blastx -p 4 -q bin.2.fa -f 6 -d uniprot_sprot.dmnd -o bin2_diamond.txt
```

!!! question
    How did that go? Did we find anything meaningful?

It is very possible the swissprot database is too small for finding meaningful hits for undersequenced / poorly known organisms.

Let us try with another piece of software, kraken2

```bash
curl -O https://ccb.jhu.edu/software/kraken2/dl/minikraken2_v1_8GB.tgz
tar xf minikraken2_v1_8GB.tgz 
```

and then

```bash
kraken2 --memory-mapping --db . --threads 4 --output bin2_kraken.txt --report bin2_kraken_report.txt metabat/bin.2.fa
```

!!! question
    Does this confirm our initial diagnostic? What kind of organism do we have?

## Functional annotation

Now that we have at least a genus for our organism, let us try to look at what it does:

```bash
# TODO prokka download
```

then we download the Actinobacteria models necessary for the functional annotation

```bash
mkdir /home/student/miniconda3/lib/python2.7/site-packages/data
download_eggnog_data.py actNOG
```

then we add functional categories to our genes:

```bash
emapper.py -d actNOG -i annotation/bin.2/bin.2.faa --no_refine -o bin2_NOG
```

and we download `bin2_NOG.emapper.annotations` to our own computers

## rpoB phylogeny

!!! warning
    We're about to produce a gene tree.
    Phylogenetic trees based on one gene must be interpreted very carefully, and may not reflect the actual phylogeny of the species

*Coming soon*


# Introduction to Nanopore Sequencing

In this tutorial we will assemble the _E. coli_ genome using a mix of long, error-prone reads from the MinION (Oxford Nanopore) and short reads from a HiSeq instrument (Illumina).

The MinION data used in this tutorial come a test run by the [Loman lab](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).  
The Illumina data were simulated using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)

## Get the Data

First download the nanopore data

```bash
wget http://s3.climb.ac.uk/nanopore/ecoli_allreads.fasta
```

You will not need the HiSeq data right away, but you can start the download in another window

```bash
curl -O -J -L https://osf.io/pxk7f/download
curl -O -J -L https://osf.io/zax3c/download
```

look at basic stats of the nanopore reads

```bash
assembly-stats ecoli_allreads.fasta
```

!!! question
How many nanopore reads do we have?

!!! question
How long is the longest read?

!!! question
What is the average read length?

## Adapter trimming

The guppy basecaller, i.e. the program that transform raw electrical signal in fastq files, already demultiplex and trim for us.

## Assembly

We assemble the reads using wtdbg2 (version > 2.3)

```bash
head -n 20000 ecoli_allreads.fasta > subset.fasta
wtdbg2 -x rs -i subset.fasta -fo assembly
wtpoa-cns -i assembly.ctg.lay -fo assembly.ctg.fa
```

## Polishing

Since the assembly likely contains a lot of errors, we correct it with Illumina reads.

First we map the short reads against the assembly

```
bowtie2-build assembly.ctg.fa assembly
bowtie2 -x assembly -1 ecoli_hiseq_R1.fastq.gz -2 ecoli_hiseq_R2.fastq.gz | \
    samtools view -bS -o assembly_short_reads.bam
samtools sort assembly_short.bam -o assembly_short_sorted.bam
samtools index assembly_short_sorted.bam
```

then we run the consensus step

```
samtools view assembly_short_sorted.bam | ./wtpoa-cns -t 16 -x sam-sr \
    -d assembly.ctg.fa -i - -fo assembly_polished.fasta
```

which will correct eventual misamatches in our assembly and write the new improved assembly to `assembly_polished.fasta`

For better results we should perform more than one round of polishing.

## Compare with the existing assembly and an illumina only assembly

### an existing assembly

Go to [https://www.ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov) and search for NC_000913.
Download the associated genome in fasta format and rename it to `ecoli_ref.fasta`

```bash
nucmer --maxmatch -c 100 -p ecoli assembly_polished.fasta ecoli_ref.fasta
mummerplot --fat --filter --png --large -p ecoli ecoli.delta
```

then take a look at `ecoli.png`

### compare metrics

!!! note
First you need to assemble the illumina data

Then run busco and quast on the 3 assemblies

!!! question
which assembly would you say is the best?

## Annotation

If you have time, train your annotation skills by running prokka on your genome!

```bash
prokka --outdir annotation --kingdom Bacteria assembly_polished.fasta
```

You can open the output to see how it went

```bash
cat annotation/*.txt
```

!!! question
Does it fit your expectations? How many genes were you expecting?

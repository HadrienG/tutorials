# Introduction to Nanopore Sequencing

In this tutorial we will assemble the *E. coli* genome using a mix of long, error-prone reads from the MinION (Oxford Nanopore) and short reads from a HiSeq instrument (Illumina).

The MinION data used in this tutorial come a test run by the [Loman lab](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).    
The Illumina data were simulated using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)

## Get the Data

First download the nanopore data

```bash
fastq-dump ERR1147227
```

You will not need the HiSeq data right away, but you can start the download in another window

```bash
curl -O -J -L https://osf.io/pxk7f/download
curl -O -J -L https://osf.io/zax3c/download
```

!!! question
    How many nanopore reads do we have?

## Adapter trimming

We'll use [porechop](https://github.com/rrwick/Porechop) to remove the adapters from the reads.
Additionally to trim the adapters at the 3' and 5' ends, porechop can split the reads if it finds adapters in the middle.

```bash
# export PATH="$PATH:/Library/Frameworks/Python.framework/Versions/3.6/bin"
porechop -i ERR1147227.fastq -o ERR1147227_trimmed.fastq
```

## Assembly

We assemble the reads using miniasm

```bash
minimap2 -x ava-ont ERR1147227_trimmed.fastq ERR1147227_trimmed.fastq | \
    gzip -1 > ERR1147227.paf.gz
miniasm -f ERR1147227_trimmed.fastq ERR1147227.paf.gz > ERR1147227.gfa
awk '/^S/{print ">"$2"\n"$3}' ERR1147227.gfa | fold > ERR1147227.fasta
```

!!! note
    Miniasm is a fast but has no consensus step.
    The accuracy of the assembly will be equal to the base accuracy.

## Polishing

Since the miniasm assembly likely contains a lot if errors, we correct it with  Illumina reads.

First we map the short reads against the assembly

```
bowtie2-build ERR1147227.fasta ERR1147227
bowtie2 -x ERR1147227 -1 ecoli_hiseq_R1.fastq.gz -2 ecoli_hiseq_R2.fastq.gz | \
    samtools view -bS -o ERR1147227.bam
samtools sort ERR1147227.bam -o ERR1147227.sorted.bam
samtools index ERR1147227.sorted.bam
```

then we run Pilon

```
pilon --genome ERR1147227.fasta --frags ERR1147227.sorted.bam \
    --output ERR1147227_improved
```

which will correct eventual misamatches in our assembly and write the new improved assembly to `ERR1147227_improved.fasta`

For better results we should perform more than one round of polishing.

## Compare with the existing assembly

Go to [https://www.ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov) and search for NC_000913.
Download the associated genome in fasta format and rename it to `ecoli_ref.fasta`

```bash
nucmer --maxmatch -c 100 -p ecoli ERR1147227_trimmed.fastq ecoli_ref.fasta
mummerplot --fat --filter --png --large -p ecoli ecoli.delta
```

then take a look at `ecoli.png`

## Annotation

```bash
awk '/^>/{print ">ctg" ++i; next}{print}' < ERR1147227_improved.fasta \
    > ERR1147227_formatted.fasta
prokka --outdir annotation --kingdom Bacteria ERR1147227_formatted.fasta
```

You can open the output to see how it went

```bash
cat annotation/PROKKA_11232017.txt
```

!!! question
    Does it fit your expecations? How many genes were you expecting?

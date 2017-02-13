# Whole-genome de-novo Assembly

In this practical we will perform the assembly of M. genitalium, a bacterium published in 1995 by Fraser et al in Science.

## Getting the data

Go to ENA, and search for the run ERR486840.

Download the 2 fastq files associated with the run.

## Quality control

Perform quality control of the data as you did in the [QC tutorial](qc.md)

How many reads are in the fastq file? What is the read length?
Does the data need trimming or other filtering? If so, do it.

Find the genome size of M. genitalium in the Fraser paper abstract.
Based on the expected genome size, the read length and the number of reads – what average coverage do you expect to get from this fastq read files?

## De-novo assembly

We will be using the SPAdes assembler to assemble our bacterium.

```
module load spades
spades.py -k 21,33,55,77,99 --careful --only-assembler -1 read_1.fastq -2 read_2.fastq -o output
```

This will produce a series of outputs. The scaffolds will be in fasta format.

How well does the assembly total consensus size and coverage correspond to your earlier estimation?
What is the N50 of the assembly? What does this mean?
How many contigs in total did the assembly produce?
How many contigs longer than 500bp? What is the N50 of those contigs only?

Perform another assembly with the following options:

Use the raw reads (no trimming, but with adapters removed), wthout the --only-assembler option.

If you have time, try out the following genome assemblers:


* MaSurCa
* Ray

## Comparing assemblies

QUAST is a software evaluating the quality of genome assemblies by computing various metrics, including

* N50, length for which the collection of all contigs of that length or longer covers at least 50% of assembly length
* NG50, where length of the reference genome is being covered
* NA50 and NGA50, where aligned blocks instead of contigs are taken
* misassemblies, misassembled and unaligned contigs or contigs bases
* genes and operons covered

To run Quast:

```
module load quast
quast.py assembly1.fasta assembly2.fasta ... -R reference.fasta -G reference.gff
```

Quast will produce a pdf in the `quast_results` directory. Download it on your computer and take a look. Which assembly is better?

## Fixing misassemblies

Pilon is a software tool which can be used to automatically improve draft assemblies. It attempts to make improvements to the input genome, including:

* Single base differences
* Small Indels
* Larger Indels or block substitution events
* Gap filling
* Identification of local misassemblies, including optional opening of new gaps

Pilon then outputs a FASTA file containing an improved representation of the genome from the read data and an optional VCF file detailing variation seen between the read data and the input genome.

You can read how Pilon works in detail [here](https://github.com/broadinstitute/pilon/wiki/Methods-of-Operation)

Before running Pilon itself, you have to map your reads back to the assembly!

```
bowtie2-build $assembly $output/index
(bowtie2 -x $output/index -1 $r1 -2 $r2 | samtools view -bS -o $mapping - ) 2> bowtie.err
samtools sort $mapping $mapping.sorted
samtools index $mapping.sorted.bam
```

Run Pilon with the following command:

```
module load pilon
pilon --genome $assembly --frags $mapping.sorted.bam --output $output
```

Once Pilon is finished running, compare the new assembly with the old one using Quast!

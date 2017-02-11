# Whole-genome de-novo Assembly

In this practical we will perform the assembly of M. genitalium, a bacterium published in 1995 by Fraser et al in Science.

## Getting the data

Go to ENA, and search for the run ERR486840.

Download the 2 fastq files associated with the run.

## Quality control?

How many reads are in the fastq file? What is the read length?
Does the data need trimming or other filtering? If so, do it.

Find the genome size of M. genitalium in the Fraser paper abstract.
Based on the expected genome size, the read length and the number of reads – what average coverage do you expect to get from this fastq read files?

## De-novo assembly

We will be using the SPAdes assembler to assemble our bacterium.

This will produce a series of outputs. The scaffolds will be in fasta format.

How well does the assembly total consensus size and coverage correspond to your earlier estimation?
What is the N50 of the assembly? What does this mean?
How many contigs in total did the assembly produce?
How many contigs longer than 500bp? What is the N50 of those contigs only?

Perform more assemblies with the following options:

Raw reads (no trimming, but with adapters removed), let spades do the qc.

## Comparing assemblies

quast

### Comparing to the reference?

Get the reference fasta here

## Fixing misassemblies?

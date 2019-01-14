# De-novo Genome Assembly

## Lecture

<iframe src="https://docs.google.com/presentation/d/e/2PACX-1vRKVI_pHGubDWeRPaAO7c9g55DzHMO5Lgd7g7AZXvjB77wAAb-wED82lXgV5P7GPF02k-21YMx8ObaX/embed?start=false&loop=false&delayms=3000" frameborder="0" width="480" height="389" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

## Practical

In this practical we will perform the assembly of _M. genitalium_, a bacterium published in 1995 by Fraser et al in Science ([abstract link](https://www.ncbi.nlm.nih.gov/pubmed/7569993)).

## Getting the data

_M. genitalium_ was sequenced using the MiSeq platform (2 \* 150bp).
The reads were deposited in the ENA Short Read Archive under the accession [ERR486840](https://www.ebi.ac.uk/ena/data/view/ERR486840)

Download the 2 fastq files associated with the run.

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_2.fastq.gz
```

The files that were deposited in ENA were already trimmed, so we do not have to trim ourselves!

!!! question
How many reads are in the files?

## De-novo assembly

We will be using the MEGAHIT assembler to assemble our bacterium

```bash
megahit -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz -o m_genitalium
```

This will take a few minutes.

The result of the assembly is in the directory m_genitalium under the name `final.contigs.fa`

Let's make a copy of it

```bash
cp m_genitalium/final.contigs.fa m_genitalium.fasta
```

and look at it

```bash
head m_genitalium.fasta
```

## Quality of the Assembly

QUAST is a software evaluating the quality of genome assemblies by computing various metrics, including

Run Quast on your assembly

```bash
quast.py m_genitalium.fasta -o m_genitalium_report
```

and take a look at the text report

```bash
cat m_genitalium_report/report.txt
```

You should see something like

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    m_genitalium
# contigs (>= 0 bp)         17
# contigs (>= 1000 bp)      8
# contigs (>= 5000 bp)      7
# contigs (>= 10000 bp)     6
# contigs (>= 25000 bp)     5
# contigs (>= 50000 bp)     2
Total length (>= 0 bp)      584267
Total length (>= 1000 bp)   580160
Total length (>= 5000 bp)   577000
Total length (>= 10000 bp)  570240
Total length (>= 25000 bp)  554043
Total length (>= 50000 bp)  446481
# contigs                   11
Largest contig              368542
Total length                582257
GC (%)                      31.71
N50                         368542
N75                         77939
L50                         1
L75                         2
# N's per 100 kbp           0.00
```

which is a summary stats about our assembly.
Additionally, the file `m_genitalium_report/report.html`

You can either download it and open it in your own web browser, or we make it available for your convenience:

- [m_genitalium_report/report.html](data/fastqc/report.html)

!!! note
N50: length for which the collection of all contigs of that length or longer covers at least 50% of assembly length

!!! question
How well does the assembly total consensus size and coverage correspond to your earlier estimation?

!!! question
How many contigs in total did the assembly produce?

!!! question
What is the N50 of the assembly? What does this mean?

## Fixing misassemblies

Pilon is a software tool which can be used to automatically improve draft assemblies.
It attempts to make improvements to the input genome, including:

- Single base differences
- Small Indels
- Larger Indels or block substitution events
- Gap filling
- Identification of local misassemblies, including optional opening of new gaps

Pilon then outputs a FASTA file containing an improved representation of the genome from the read data and an optional VCF file detailing variation seen between the read data and the input genome.

Before running Pilon itself, we have to align our reads against the assembly

```
bowtie2-build m_genitalium.fasta m_genitalium
bowtie2 -x m_genitalium -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz | \
    samtools view -bS -o m_genitalium.bam
samtools sort m_genitalium.bam -o m_genitalium.sorted.bam
samtools index m_genitalium.sorted.bam
```

then we run Pilon

```
pilon --genome m_genitalium.fasta --frags m_genitalium.sorted.bam --output m_genitalium_improved
```

which will correct eventual mismatches in our assembly and write the new improved assembly to `m_genitalium_improved.fasta`

## Assembly Completeness

Although quast output a range of metric to assess how contiguous our assembly is, having a long N50 does not guarantee a good assembly: it could be riddled by misassemblies!

We will run `busco` to try to find marker genes in our assembly. Marker genes are conserved across a range of species and finding intact conserved genes in our assembly would be a good indication of its quality

First we need to download and unpack the bacterial datasets used by `busco`

```bash
wget http://busco.ezlab.org/datasets/bacteria_odb9.tar.gz
tar xzf bacteria_odb9.tar.gz
```

then we can run `busco` with

```bash
BUSCO.py -i m_genitalium.fasta -l bacteria_odb9 -o busco_genitalium -m genome
```

!!! question
How many marker genes has `busco` found?

## Course literature

Course litteraturer for today is:

- Next-Generation Sequence Assembly: Four Stages of Data Processing and Computational Challenges: <https://doi.org/10.1371/journal.pcbi.1003345>

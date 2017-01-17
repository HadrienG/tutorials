# File Formats

This lecture is aimed at making you discover the most popular file formats used in bioinforatics. You're expected to have basic working knowledge of linux to be able to follow the lesson.

## Table of Contents
* [The fasta format](#the-fasta-format)
* [The fastq format](#the-fastq-format)
* [The sam/bam format](#the-sam-format)
* [The vcf format](#the-vcf-format)
* [The gff format](#the-gff-format)

### The fasta format

The fasta format was invented in 1988 and designed to represent nucleotide or peptide sequences. It originates from the [FASTA](https://en.wikipedia.org/wiki/FASTA) software package, but is now a standard in the world of bioinformatics

The first line in a FASTA file starts with a ">" (greater-than) symbol followed by the description or identifier of the sequence. Following the initial line (used for a unique description of the sequence) is the actual sequence itself in standard one-letter code.

A few sample sequences:

```
>KX580312.1 Homo sapiens truncated breast cancer 1 (BRCA1) gene, exon 15 and partial cds
GTCATCCCCTTCTAAATGCCCATCATTAGATGATAGGTGGTACATGCACAGTTGCTCTGGGAGTCTTCAG
AATAGAAACTACCCATCTCAAGAGGAGCTCATTAAGGTTGTTGATGTGGAGGAGTAACAGCTGGAAGAGT
CTGGGCCACACGATTTGACGGAAACATCTTACTTGCCAAGGCAAGATCTAG
```

```
>KRN06561.1 heat shock [Lactobacillus sucicola DSM 21376 = JCM 15457]
MSLVMANELTNRFNNWMKQDDFFGNLGRSFFDLDNSVNRALKTDVKETDKAYEVRIDVPGIDKKDITVDY
HDGVLSVNAKRDSFNDESDSEGNVIASERSYGRFARQYSLPNVDESGIKAKCEDGVLKLTLPKLAEEKIN
GNHIEIE
```

A fasta file can contain multiple sequence. Each sequence will be separated by their "header" line, starting by ">"

Example:

```
>KRN06561.1 heat shock [Lactobacillus sucicola DSM 21376 = JCM 15457]
MSLVMANELTNRFNNWMKQDDFFGNLGRSFFDLDNSVNRALKTDVKETDKAYEVRIDVPGIDKKDITVDY
HDGVLSVNAKRDSFNDESDSEGNVIASERSYGRFARQYSLPNVDESGIKAKCEDGVLKLTLPKLAEEKIN
GNHIEIE
>3HHU_A Chain A, Human Heat-Shock Protein 90 (Hsp90)
MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNKEIFLRELISNSSDALDKIRYESLTDPSKL
DSGKELHINLIPNKQDRTLTIVDTGIGMTKADLINNLGTIAKSGTKAFMEALQAGADISMIGQFGVGFYS
AYLVAEKVTVITKHNDDEQYAWESSAGGSFTVRTDTGEPMGRGTKVILHLKEDQTEYLEERRIKEIVKKH
SQFIGYPITLFVEK
```

### The fastq format

The fastq format is also a text based format to represent nucleotide sequences, but also contains the coresponding quality of each nucleotide. It is the standard for storing the output of high-throughput sequencing instruments such as the Illumina machines.

A fastq file uses four lines per sequence:

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

An example sequence in fastq format:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

#### Quality

The quality, also called phred scores is the probability that the corresponding basecall is incorrect.

Phred scores use a logarithmic scale, and are represented by ASCII characters, mapping to a quality usually going from 0 to 40

Phred Quality Score | Probability of incorrect base call | Base call accuracy
- | - | -
10 | 1 in 10 | 90%
20 | 1 in 100 | 99%
30 | 1 in 1000 | 99.9%
40 | 1 in 10,000 | 99.99%
50 | 1 in 100,000 | 99.999%
60 | 1 in 1,000,000 | 99.9999%

### the sam/bam format



### the vcf format

### the gff format

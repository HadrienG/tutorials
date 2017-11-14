# Quality control

In this practical you will learn to import, view and check the quality of NGS read data in FASTQ format.

You will be working with an Illumina MiSeq read dataset from a genome sequence project. The sequenced organism is an enterohaemorrhagic E. coli (EHEC) of the serotype O157, a potentially fatal gastrointestinal pathogen. The sequenced bacterium was part of an outbreak investigation in the St. Louis area, USA in 2011.
The sequencing was done as paired-end 2x150bp.

## Downloading the data

The raw data were deposited at the European Nucleotide Archive, under the accession number SRR957824. Go to the ENA [website](http://www.ebi.ac.uk/ena) and search for the run with the accession SRR957824. Download the two fastq files associated with the run:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR957/SRR957824/SRR957824_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR957/SRR957824/SRR957824_2.fastq.gz
```

## FastQC

To check the quality of the sequence data we will use a tool called FastQC. With this you can check things like read length distribution, quality distribution across the read length, sequencing artifacts and much more.

FastQC has a graphical interface and can be downloaded and run on a Windows or Linux computer without installation. It is available [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

However, FastQC is also available as a command line utility on the training server you are using. You can load the module and execute the program as follows:

```
module load FastQC
fastqc $read1 $read2
```

which will produce both a .zip archive containing all the plots, and a html document for you to look at the result in your browser.

Open the html file with your favourite web browser, and try to interpret them. 

Pay special attention to the per base sequence quality and sequence length distribution. Explanations for the various quality modules can be found [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/). Also, have a look at examples of a [good](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and a [bad](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) illumina read set for comparison.

You will note that the reads in your uploaded dataset have fairly poor quality (<20) towards the end. There are also outlier reads that have very poor quality for most of the second half of the reads.

There are overrepresented sequences in the data. Where do they come from?

## Scythe

Scythe uses a Naive Bayesian approach to classify contaminant substrings in sequence reads. It considers quality information, which can make it robust in picking out 3'-end adapters, which often include poor quality bases.

First, install scythe:

```
git clone https://github.com/vsbuffalo/scythe.git
cd scythe
make all
```

Then, copy or move "scythe" to a directory in your $PATH, for example like this:

`cp scythe $HOME/bin/`

Scythe can be run minimally with:

`scythe -a adapter_file.fasta -o trimmed_sequences.fastq sequences.fastq`

Try to trim the adapters in both your read files!

## Sickle

Most modern sequencing technologies produce reads that have deteriorating quality towards the 3'-end and some towards the 5'-end as well. Incorrectly called bases in both regions negatively impact assembles, mapping, and downstream bioinformatics analyses.

We will trim each read individually down to the good quality part to keep the bad part from interfering with downstream applications.

To do so, we will use sickle. Sickle is a tool that uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim the 3'-end of reads and also determines when the quality is sufficiently high enough to trim the 5'-end of reads. It will also discard reads based upon a length threshold.

First, install sickle:

```
git clone https://github.com/najoshi/sickle.git
cd sickle
make
```

Copy sickle to a directory in your $PATH:

`cp sickle $HOME/bin/`

Sickle has two modes to work with both paired-end and single-end reads: sickle se and sickle pe.

Running sickle by itself will print the help:

`sickle`

Running sickle with either the "se" or "pe" commands will give help specific to those commands. Since we have paired end reads:

`sickle pe`

Set the quality score to 25. This means the trimmer will work its way from both ends of each read, cutting away any bases with a quality score < 25.

```
sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
-o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
-s trimmed_singles_file.fastq -q 25
```

What did the trimming do to the per-base sequence quality, the per sequence quality scores and the sequence length distribution? Run FastQC again to find out.

What is the sequence duplication levels graph about? Why should you care about a high level of duplication, and why is the level of duplication very low for this data?

Based on the FastQC report, there seems to be a population of shorter reads that are technical artifacts. We will ignore them for now as they will not interfere with our analysis.

## Extra exercises

Perform quality control on the extra datasets given by your instructors. 

# Mapping and variant calling

In this practical you will learn to map NGS reads to a reference sequence, check the output using a viewer software and investigate some aspects of the results. You will be using the read data from the [Quality Control](qc.md) practical.

EHEC O157 strains generally carry a large virulence plasmid, pO157. Plasmids are circular genetic elements that many bacteria carry in addition to their chromosomes. This particular plasmid encodes a number of proteins which are known or suspected to be involved in the ability to cause severe disease in infected humans. Your task in this practical is to map your prepared read set to a reference sequence of the virulence plasmid, to determine if the pO157 plasmid is present in the St. Louis outbreak strain.

## Read mapping

### downloading a reference

You will need a reference sequence to map your reads to. You can download the reference from [here](data/pO157_Sakai.fasta)

or from the command-line with

`wget https://github.com/HadrienG/tutorials/tree/master/data/pO157_Sakai.fasta`

This file contains the sequence of the pO157 plasmid from the Sakai outbreak strain of E. coli O157. In contrast to the strain we are working on, this strain is available as a finished genome, i.e. the whole sequence of both the single chromosome and the large virulence plasmid are known.

### indexing the reference

```
module load bowtie2
bowtie2-build $ref $prefix
```

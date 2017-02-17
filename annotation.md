# Genome Annotation

After you have de novo assembled your genome sequencing reads into contigs, it is useful to know what genomic features are on those contigs. The process of identifying and labelling those features is called genome annotation.

Prokka is a “wrapper”; it collects together several pieces of software (from various authors), and so avoids “re-inventing the wheel”.

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence. Prokka uses a two-step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://prodigal.ornl.gov); second, the function of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

## Input data

Prokka requires assembled contigs. You will need your best assembly from the assembly tutorial.

Alternatively, you can download an assembly [here](data/assembly.fasta)

## Running prokka

```
module load prokka
awk '/^>/{print ">ctg" ++i; next}{print}' < assembly.fasta > good_contigs.fasta
prokka --outdir annotation --kingdom Bacteria \
--proteins m_genitalium.faa good_contigs.fasta
```

Once Prokka has finished, examine each of its output files.

* The GFF and GBK files contain all of the information about the features annotated (in different formats.)
* The .txt file contains a summary of the number of features annotated.
* The .faa file contains the protein sequences of the genes annotated.
* The .ffn file contains the nucleotide sequences of the genes annotated.

## Visualising the annotation

Artemis is a graphical Java program to browse annotated genomes. Download it [here](http://www.sanger.ac.uk/science/tools/artemis) and install it on your local computer.

Copy the .gff file produced by prokka on your computer, and open it with artemis.

You will be overwhelmed and/or confused at first, and possibly permanently. Here are some tips:

* There are 3 panels: feature map (top), sequence (middle), feature list (bottom)
* Click right-mouse-button on bottom panel and select Show products
* Zooming is done via the verrtical scroll bars in the two top panels

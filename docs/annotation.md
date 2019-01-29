# Genome Annotation

## Lecture

<br>

<iframe src="https://docs.google.com/presentation/d/e/2PACX-1vTERGc6gJyJeGylr6xzXvioMFixfI6x9XIT8QHqC8XIq8cP3KHe6PUuumbMrunSCVlbFhFJaVh2wvMh/embed?start=false&loop=false&delayms=3000" frameborder="0" width="480" height="389" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

After you have de novo assembled your genome sequencing reads into contigs, it is useful to know what genomic features are on those contigs. The process of identifying and labelling those features is called genome annotation.

Prokka is a "wrapper"; it collects together several pieces of software (from various authors), and so avoids "re-inventing the wheel".

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence. Prokka uses a two-step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://compbio.ornl.gov/prodigal/); second, the function of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

## Input data

Prokka requires assembled contigs. You can prepare you working directory for this annotation tutorial.

```bash
mkdir ~/annotation
cd ~/annotation
```

You will download an improved assembly of *Mycoplasma genitalium* into you data directory:

```bash
curl -O -J -L https://osf.io/7eaky/download
```

You will also need a proteins set specific of Mycoplasma for the annotation. Here is a file containing the Mycoplasma proteins retrieved from Swiss-Prot database (3041 sequences)

```bash
curl -O -J -L https://osf.io/xjm3n/download
```


## Running prokka

```bash
prokka --outdir annotation --kingdom Bacteria \
--proteins uniprot_mycoplasma_reviewed.faa m_genetalium_improved.fasta
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

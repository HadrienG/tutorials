# Pan-Genome Analysis

In this tutorial we will learn how to determine a pan-genome from a collection of isolate genomes.

This tutorial is inspired from [Genome annotation and Pangenome Analysis](https://github.com/microgenomics/tutorials/blob/master/pangenome.md) from the CBIB in Santiago, Chile

## Getting the data

We'll data from [this article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006258) and analyse the core and accessory genomes of _E. coli strains_

Firstly, download the supplementary csv file containing information on all the strains using in the study.

- [Link](https://doi.org/10.1371/journal.pcbi.1006258.s010)

Then, open it in Excel (or any software that opens .csv files) and select 32 strains.
Download these 32 strains from the ENA website!

!!! note
    Alternatively you can use [ena browser tools](https://github.com/enasequence/enaBrowserTools) to download the files.
    It is available as a bioconda recipe.

```bash
# for getting 10 random strains at the command-line
cut -d',' -f 1 journal.pcbi.1006258.s010.csv | tail -n +2 | shuf | head -10 > strains.txt
cat strains.txt | parallel enaGroupGet -f fastq {}
```

if you wish to use the same strains as your instructor:


```bash
curl -O -J -L https://osf.io/s43mv/download
cat strains.txt | parallel enaGroupGet -f fastq {}
```

and then put all the reads in the same directory

```bash
mkdir reads
mv ERS*/*/*.fastq.gz reads/
rm -r ERS*
```

## Assemble and Annotate the strains

You'll assemble your strains with megahit.

```bash
mkdir assemblies
for r1 in reads/*_1.fastq.gz
do
    prefix=$(basename $r1 _1.fastq.gz)
    r2=reads/${prefix}_2.fastq.gz
    megahit -1 $r1 -2 $r2 -o ${prefix} --out-prefix ${prefix}
    mv ${prefix}/${prefix}.contigs.fa assemblies/
    rm -r ${prefix}
done
```

and use prokka to annotate


```bash
mkdir annotation
for assembly in assemblies/*.fa
do
    prefix=$(basename $assembly .contigs.fa)
    prokka --usegenus --genus Escherichia --species coli --strain ${prefix} \
        --outdir ${prefix} --prefix ${prefix} ${assembly}
    mv ${prefix}/${prefix}.gff annotation/
    rm -r ${prefix}
done
```

## Pan-genome analysis

put all the .gff files in the same folder (e.g., `./gff`) and run Roary

```bash
roary -f roary -e -n -v annotation/*.gff
```

Roary will get all the coding sequences, convert them into protein, and create pre-clusters. Then, using BLASTP and MCL, Roary will create clusters, and check for paralogs. Finally, Roary will take every isolate and order them by presence/absence of orthologs. The summary output is present in the `summary_statistics.txt` file.

Additionally, Roary produces a `gene_presence_absence.csv` file that can be opened in any spreadsheet software to manually explore the results. In this file, you will find information such as gene name and gene annotation, and, of course, whether a gene is present in a genome or not.

## Plotting the result

Roary comes with a python script that allows you to generate a few plots to graphically assess your analysis output.

First, we need to generate a tree file from the alignment generated by Roary:

```
cd roary
FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick
```

Then we can plot the Roary results with `roary_plots.py`, a community contriubuted python script to visualise roary results:

```
wget https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py
python roary_plots.py
roary_plots.py my_tree.newick gene_presence_absence.csv
```

then look at the 3 `/png` files that have been generated

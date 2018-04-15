# Structuring a python project

In this tutorial you will learn how to structure a python project.
In order to follow this you will need basic knowledge of git.

We'll create a small python package step by step.
Our package will give us statistics about genome sequences or reads, as well as convert common file formats used in bioinformatics.

!!! note
    `nano` will start to show its limits as a file editor. We recommend the excellent `atom` for this lesson.
    Additionally it has built-in git support!

!!! note
    since this course assumes you have git knowledge, you'll be expected to commit everytime you see a not "time to commit!"


## Creating your project

We will call you project `mfpp`, short for my first python package!

```bash
mkdir mfpp
cd mffp
git init
```

## The structure of a project

Inside your newly created mfpp directory, our finished python project will look like this:

```
.
..
.git
mfpp/
    __init__.py
    __main__.py
    app.py
    util.py
    reader.py
    test/
        test_util.py
        test_reader.py
requirements.txt
setup.py
LICENSE
README.md
```

Don't worry, we'll explain what each file does, and create them along the way

## Modules

We want python to see our project as a module, meaning that we can execute our program at the command line using `python -m mfpp`

For that we need a few things. Let us first create (again) a directory called `mfpp`

### \_\_init\_\_.py

```python
from mfpp import *
```

### \_\_main\_\_.py

```python
from mfpp.app import main
main()
```

### app.py

This is the core of our project. it should contain a `main` function

```python
def main():
    print('hello world!')
```

and now in the terminal:

```bash
python -m mfpp
# hello world!
```

the only thing we need to have a fully fledged python package (besides a license, documenation, ...) is `setup.py`
setup.py will allow our package to be installable with `pip`, by anyone on any computer equipped with python3
We will not create setup.py right away, since we can test the functionalities of our package with `python -m mpff`

!!! Reminder
    Time to commit!

## Our first functionality

in this chapter we will add a functionaily to our `mfpp` reading a fasta file and output some stats about it

the finished product:

```bash
python -m mfpp --input my_file.fasta
```

### Libraries

for some parts of our program, such as reading a fasta file, some people have already created python libraries that do what we want to do.

To not re-invent the wheel, we will use those libraries.
But we do not want them to install them globally on our computer, we'll use a virtual environment instead

```
pip3 install virtualenv
virtualenv venv_mfpp
source venv_mfpp/bin/activate
```

!!! note
    use deactivate to exit the virtual environment

so we do not pollute our system python with libraries

we do not want to keep track of our virtualenv in git so we create a `.gitignore` and add it to it.

```
# gitignore
*.pyc
venv_mfpp
```

Now we can install the needed libraries

```
# requirements.txt
biopython
```

```bash
pip install -r requirements.txt
```

!!! Reminder
    Time to commit!

### reading a fasta file with biopython

First download the *E. coli strain k12* genome and put in in a data directory

```bash
mkdir data
wget $todo_path_to_genome
```

then in app.py

```python
from Bio import SeqIO


def main():
    genome_path = 'data/ecoli.fasta'
    with open(genome_path, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            print(record.id)
```

```bash
python -m mfpp
# NC_000913.3
```

it works! but it always work on the same fast file, which is not very interesting.

!!! Reminder
    Time to commit!

### Command line arguments with `argparse`

```python
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        prog='mfpp',
        description='my first python project'
    )
    parser.add_argument(
        '--input',
        type=str,
        help='an input fasta file'
    )
    args = parser.parse_args()

    genome_path = args.input
    with open(genome_path, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            print(record.id)
```

we now have a help message!

```bash
python -m mpff -h
usage: mfpp [-h] input

my first python project

positional arguments:
  input       an input fasta file

optional arguments:
  -h, --help  show this help message and exit
```

and if we execute

```bash
python -m mfpp data/ecoli.fasta
# NC_000913.3
```

As our application grows, we may have more and more arguments in our main.
Therefore it is good practise to move the functionailties of our software in other files, and to import them.

Let us move our fasta parsin code in a new file called `fasta_stats.py` and reexecute it.

## Writing tests

## The README

## A word on Licenses

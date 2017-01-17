# Introduction to Unix (continued)

In this part of the tutorial, we'll learn how to install programs in a Unix system

## Using a package manager

This is the most straight-forward way, and the way used by most of the people using unix at home,
or administrating their own machine.

This course is aimed at giving you a working knowledge of linux for bioinformatics, and in that setting, you will rarely, if ever, be the administrator of your own machine. The methods below are here as an information

### On Ubuntu and Debian: Apt

To install a software:

`apt-get install name_of_the_software`

to uninstall:

`apt-get remove name_of_the_software`

to update all installed softwares:

```
apt-get update
apt-get upgrade
```

### On Fedora, CentOS and RedHat: yum

To install a software:

`yum install name_of_the_software`

to uninstall:

`yum remove name_of_the_software`

to update:

`yum update`

### MacOS: brew

Although there are no official package managers on MacOS, two popular, community-driven alternatives exist: macports and brew.

Brew is particularly pupular within the bioinformatics community, and allows easy installation of many bioinformatics softwares on MacOS

To install brew on your mac:

`/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`

To install a software:

`brew install name_of_the_software`

To uninstall:

`brew uninstall name_of_the_software`

To update all brew-installed softwares:

```
brew update
brew upgrade
```

More info on [brew.sh](http://brew.sh) and [brew.sh/homebrew-science/](http://brew.sh/homebrew-science/)

## Downloading binaries

In a university setting, you will rarely by administrator of your own machine. This is a very good thing for one reason: it's harder for you to break something!

The downside is that it makes installing softwares more complicated. We'll start wit simply downloading the software and executing it, then we'll learn how to obtain packages from source code.

for example, we'll install the blast binaries:

First, download the archive:
`wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz`

then unpack it and go to the newly created directory

```
tar xzf ncbi-blast-2.6.0+-x64-linux.tar.gz
cd ncbi-blast-2.6.0+
```

you should have a `bin` directory, go inside and look at the files. You have a bunch of executable files.

### Execute a file


Most of the lunix commands that you execute on a regular basis (ls, cp, mkdir) are located in `/usr/bin`, but you don't have to invoke them with their full path: i.e. you dont type `/usr/bin/ls` but just `ls`. This is because `/usr/bin/ls` is in your $PATH.

to execute a file that you just downloaded, and is therefore not in your path, you have to type the absolute or relative path to that file. Meaning, for the blast program suite that we just downloaded:

`bin/blastn -help`

or

```
cd bin
./blastn -help
```

and that's it!
But it is not very convenient. You want to be able to execute blast without having to remember where it is. If you have administrator rights (sudo), you can move the software in `/usr/bin`. If you don't you can modify your $PATH in a configuration file called `.bash_profile` that is located in your home.

More information on how to correctly modify your PATH [here](http://unix.stackexchange.com/a/26059)

## Compiling from source

Sometimes pre-compiled binaries are not available. You then have to compile from source: transforming the human-readable code (written in one or another programming language) into machine-readable code (binary)

The most common way to do so, if a software package has its source coud available online is

```
./configure
make
make install
```

If you don't have the administrator rights, you'll often have to pass an extra argument to ./configure:

```
./configure --prefix="/where_i_want_to_install"
make
make install
```

Most of the softwares come with instructions on how to install them. Always read the file called README or INSTALL in the package directory before installing!

### Exercice

The most popular unix distributions come with a version of python (a programming language) that is not the most recent one. Install from source the most recent version of python in a folder called `bin` in your home directory.

You can download the python source code at https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz

## Install python packages

Python is a really popular programming language in the world of bioinformatics. Python has a package manager called `pip` that you can use to install softwares written in python.

Please us the python executable you installed in the above exercise!

Firstly, get pip:

`wget https://bootstrap.pypa.io/get-pip.py`

then execute the script

`python get-pip.py`

Thenm you can use pip to install package, either globally (if you're an administrator):

`pip install youtube_dl`

or just for you:

`pip install --user youtube_dl`

## Final exercise

One of the oldest and most famous bioinformatics package is called EMBOSS.
Install EMBOSS in the bin directory of your home. Good luck!

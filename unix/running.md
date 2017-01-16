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

# Unix

This tutorial is largely inspired of the [Introduction to UNIX](ftp://ftp.sanger.ac.uk/pub/project/pathogens/jm15/unix.pdf) course from the Sanger Institute.

The aim of this module is to introduce Unix and cover some of the basics that will allow you to be more comfortable with the command-line. Several of the programs that you are going to use during this course are useful for bioinformatics analyses. This module is only designed to provide a very brief introduction to some of the features and useful commands of Unix. During this module we will also obtain a genome sequence and examine the basic structure of an EMBL entry.

## Introduction

Unix is the standard operating system on most large computer systems in scientific research, in the same way that Microsoft Windows is the dominant operating system on desktop PCs. Unix and MS Windows both perform the important job of managing the computer’s hardware (screen, keyboard, mouse, hard disks, network connections, etc...) on your behalf. They also provide you with tools to manage your files and to run application software. They both offer a graphical user interface (desktop). The desktops look different, call things by different names but they mostly can do the same things. Unix is a powerful, secure, robust and stable operating system that allows dozens of people to run programs on the same computer at the same time. This is why it is the preferred operating system for large-scale scientific computing. It is run on all kind of machines, like mobile phones (Android), desktop PCs, kitchen appliances,... all the way up to supercomputers. Unix powers the majority of the Internet.

## Aims

The aim of this course is to introduce Unix and cover the basics. The programs that you are going to use during the courses, plus many others that are useful for bioinformatics analyses, are run in Unix. This module is only designed to provide a very brief introduction to some of the features and useful commands of Unix. During this module we will also obtain a genome sequence and examine the basic structure of an EMBL entry.

## Why use Unix?

* Unix is a well established, very widespread operating system. You probably have a device running on Unix in your home without realising it (e.g. playstation, TV box, wireless router, android tablets/phones,...
* Command line driven, with a huge number of often terse, but powerful commands.
* In contrast to Windows, it is designed to allow many users to run their programs simultaneously on the same computer.
* Designed to work in computer networks - for example, most of the Internet is Unix based.
* It is used on many of the powerful computers at bioinformatics centres and also on many desktops and laptops (MacOS is largely UNIX compatible).
* The major difference between Unix and Windows is that it is free (as in freedom) and you can modify it to work however you want. This same principle of freedom is also used in most bioinformatics software.
* There are many distributions of Unix such as Ubuntu, RedHat, Fedora, Mint,...). These are all Unix, but they bundle up extra software in a different way or combinations. Some are known for being conservative and reliable; whilst others are know for being cutting-edge (and less reliable).
* The MacOSX operating system used by the [eBioKit](77.235.253.122) is also based on Unix.

## Getting started

For this course, you will have to connect to the eBiokit using SSH. SSH stands for Secure Shell and is a network protocol used to securely connect to a server. To do so, you will need an SSH client:

* On Linux: it is included by default, named Terminal.
* On MacOS: it is included by default, also named Terminal.
* On Windows: you'll have to download and install [MobaXterm](http://mobaxterm.mobatek.net), a terminal emulator.

Once you've opened your terminal (or terminal emulator), type

`ssh username@ip_address`

replacing `username` and `ip_address` with your username and the ip address of the server you are connecting to.
Type your password when prompted. As you type, nothing will show on screen. No stars, no dots.
It is supposed to be that way. Just type the password and press enter!

![Unix Prompt](http://77.235.253.122/tutorials/wp-content/uploads/2016/09/unix_prompt.png)

You can type commands directly into the terminal at the ‘$' prompt. A list of useful commands can be found on the next page. Many of them are two- or three-letter abbreviations. The earliest Unix systems (*circa* 1970) only had slow Teletype terminals, so it was faster to type 'rm' to remove a file than 'delete' or 'erase'. This terseness is a feature of Unix that still survives.

## The command line

All Unix programs may be run by typing commands at the Unix prompt. The command line tells the computer what to do. You may subtly alter these commands by specifying certain options when typing in the command line.

### Command line Arguments

Typing any Unix command for example **ls**, **mv** or **cd** at the Unix prompt with the appropriate variables such as files names or directories will result in the tasks being performed on pressing the enter key.

![Command Arguments](http://77.235.253.122/tutorials/wp-content/uploads/2016/09/unix_prompt_arguments.png)

The command is separated from the options and arguments by a space. Additional options and/or arguments can be added to the commands to affect the way the command works. Options usually have one dash and a letter (e.g. -h) or two dashes and a word (--help) with no space between the dash and the letter/word. Arguments are usually filenames or directories.

For example: List the contents of a directory

**ls** List the contents of a directoryList the contents of a directory with extra information about the files
**ls –l** List the contents of a directory with extra information about the files
**ls –a** List all contents including hidden files & directories
**ls -al** List all contents including hidden files & directories, with extra information about the files
**ls –l /usr/** List the contents of the directory /usr/, with extra information about the files

### Files and Directories

Directories are the Unix equivalent of folders on a PC or Mac. They are organised in a hierarchy, so directories can have sub-directories. Directories are very useful for organising your work and keeping your account tidy - for example, if you have more than one project, you can organise the files for each project into different directories to keep them separate. You can think of directories as rooms in a house. You can only be in one room (directory) at a time. When you are in a room you can see everything in that room easily. To see things in other rooms, you have to go to the appropriate door and crane your head around. Unix works in a similar manner, moving from directory to directory to access files. The location or directory that you are in is referred to as the current working directory.

**Directory structure example**

![filesystem structure](http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_3_directory.png)

Therefore if there is a file called genome.seq in the **dna** directory its location or full pathname can be expressed as /nfs/dna/genome.seq.

### General Points

Unix is pretty straightforward, but there are some general points to remember that will make your life easier: most flavors of UNIX are case sensitive - typing **ls** is generally not the same as typing **LS**. You need to put a space between a command and its argument - for example, **less my_file** will show you the contents of the file called my_file; **lessmyfile** will just give you an error! Unix is not psychic: If you misspell the name of a command or the name of a file, it will not understand you. Many of the commands are only a few letters long; this can be confusing until you start to think logically about why those letters were chosen - ls for list, rm for remove and so on. Often when you have problems with Unix, it is due to a spelling mistake, or perhaps you have omitted a space. If you want to know more about Unix and its commands there are plenty of resources available that provide a more comprehensive guide (including a cheat sheet at the end of this chapter.

*     [http://unix.t-a-y-l-o-r.com/](http://unix.t-a-y-l-o-r.com/)

In what follows, we shall use the following typographical conventions: Characters written in **bold typewriter font** are commands to be typed into the computer as they stand. Characters written in *italic typewriter font* indicate non-specific file or directory names. Words inserted within square brackets [Ctrl] indicate keys to be pressed. So, for example,
`$ **ls** *any_directory* [Enter]` means "at the Unix prompt $, type ls followed by the name of some directory, then press Enter"
Don't forget to press the [Enter] key: commands are not sent to the computer until this is done.

### Some useful Unix commands Command and What it does

| Command | What it does
--- | ---
| **ls** | Lists the contents of the current directory
| **mkdir** | Creates a new directory
| **mv** | Moves or renames a file
| **cp** | Copies a file
| **rm** | Removes a file
| **cat** | Concatenates files
| **less** | Displays the contents of a file one page at a time
| **head** | Displays the first ten lines of a file
| **tail** | Displays the last ten lines of a file
| **cd** | Changes current working directory
| **pwd** | Prints working directory
| **find** | Finds files matching an expression
| **grep** | Searches a file for patterns
| **wc** | Counts the lines, words, characters, and bytes in a file
| **kill** | Stops a process
| **jobs** | Lists the processes that are running  

### Firts steps
The following exercise introduces a few useful Unix commands and provides examples of how they can be used. Many people panic when they are confronted with an Unix prompt! Don’t! The exercise is designed to be step-by-step, so all the commands you need are provided in the text. If you get lost ask a demonstrator. If you are a person skilled at Unix, be patient it is only a short exercise. Finding where you are and what you’ve got

**pwd**     
Print the working directory As seen previously directories are arranged in a hierarchical structure. To determine where you are in the hierarchy you can use the pwd command to display the name of the current working directory. The current working directory may be thought of as the directory you are in, i.e. your current position in the file-system tree To find out where you are type

`pwd [enter]`

You will see that you are in your home directory. We need to move into the ngs_course_data directory. Remember, Unix is case sensitive **PWD** is not the same as **pwd**

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise1.png"><img class="alignnone size-full wp-image-432" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise1.png" alt="Unix_exercise1" width="865" height="177" /></a>
</p>

**cd** 
Change current working directory The cd command will change the current working directory to another, in other words allow you to move up or down in the directory hierarchy. First of all we are going to move into the "ngs_course_data" directory below. To do this type:

`cd ngs_course_data [enter]`

Now use the pwd command to check your location in the directory hierarchy.

Change again the directory to `Module_Unix`  

**ls**  
List the contents of a directory To find out what are the contents of the current directory type **ls** [enter] The ls command lists the contents of your current directory, this includes files and directories You should see that there are several other directories.

Now use the `cd` command again to change to the `Module_Unix` directory.

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise2.png"><img class="alignnone size-full wp-image-433" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise2.png" alt="Unix_exercise2" width="865" height="336" /></a>
</p>


### Changing and moving what you’ve got

**cp**  
Copy a file. **cp file1 file2** is the command which makes a copy of file1 in the current working directory and calls it file2! What you are going to do is make a copy of AL513382.embl. This file contains the genome of Salmonella typhi strain CT18 in EMBL format (we'll learn more about file formats later during the course). The new file will be called S_typhi.embl.  
`cp AL513382.embl S_typhi.embl [enter]`
If you use the ls command to check the contents of the current directory you will see that there is an extra file called S_typhi.embl.

<p style="text-align: center;">
   <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise3.png"><img class="alignnone size-full wp-image-434" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise3.png" alt="Unix_exercise3" width="865" height="206" /></a>
</p>

**rm**
Delete a file. This command removes a file permanently, so be careful! You are now going to remove the old version of the *S. typhi* genome file, AL513382.embl
`rm AL513382.embl [enter]`

The file will be removed. Use the **ls** command to check the contents of the current directory to see that AL513382.embl has been removed.

Unix, as a general rule does exactly what you ask, and does not ask for confirmation. Unfortunately there is no "recycle bin" on the command line to recover the file from, so you have to be careful.

<p style="text-align: center;">
   <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise4.png"><img class="alignnone size-full wp-image-435" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise4.png" alt="Unix_exercise4" width="865" height="185" /></a>
</p>

**cd**
Change current working directory. As before the cd command will change the current working directory to another, in other words allow you to move up or down in the directory hierarchy. First of all we are going to move into the directory above, type:
`cd .. [enter]`

Now use the **pwd** command to check your location in the directory hierarchy. Next, we are going to move into the Module_Artemis directory. To change to the Module_Artemis directory type:
`cd Module_Artemis [enter]` use the ls command to check the contents of the directory.

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise5.png"><img class="alignnone size-full wp-image-436" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise5.png" alt="Unix_exercise5" width="865" height="231" /></a>
</p>

#### Tips

There are some short cuts for referring to directories:

```
.  Current directory (one full stop)  
.. Directory above (two full stops)  
~  Home directory (tilde)  
/  Root of the file system (like C:\ in Windows)
```

Pressing the tab key twice will try and autocomplete what you’ve started typing or give you a list of all possible completions. This saves a lot of typing and typos. Pressing the up/down arrows will let you scroll through the previous commands. If you highlight some text, middle clicking will paste it on the command line.


**mv**
Move a file. To move a file from one place to another use the **mv** command. This moves the file rather than copies it, therefore you end up with only one file rather than two. When using the command the path or pathname is used to tell Unix where to find the file. You refer to files in other directories by using the list of hierarchical names separated by slashes. For example, the file bases in the directory genome has the path genome/bases If no path is specified Unix assumes that the file is in the current working directory. What you are going to do is move the file S_typhi.embl from the Module_Unix directory, to the current working directory.
`mv ../Module_Unix/S_typhi.embl . [enter]`  
Use the ls command to check the contents of the current directory to see that S_typhi.embl has been moved. ../Module_Unix/S_typhi.embl specifies that S_typhi.embl is in the Module_Unix directory. If the file was in the directory above, the path would change to: ../ S_typhi.embl

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise6.png"><img class="alignnone size-full wp-image-437" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise6.png" alt="Unix_exercise6" width="865" height="222" /></a>
</p>

The command can also be used to rename a file in the current working directory. Previously we used the cp command, but mv provides an alternative without the need to delete the original file. Therefore we could have used:  

`mv AL513382.embl S_typhi.embl [enter]` instead of:

```
cp AL513382.embl S_typhi.embl [enter]
rm AL513382.embl [enter]
```

### Viewing what you’ve got

**less**
Display file contents. This command displays the contents of a specified file one screen at a time. You are now going to look at the contents of S_typhi.embl.
`less S_typhi.embl [enter]`  
The contents of S_typhi.embl will be displayed one screen at a time, to view the next screen press the space bar. less can also scroll backwards if you hit the b key. Another useful feature is the slash key, /, to search for a word in the file. You type the word you are looking for and press enter. The screen will jump to the next occurrence and highlight it. As S_typhi.embl is a large file this will take a while, therefore you may want to escape or exit from this command. To exit press the letter ‘q’. If you really need to exit from a program and it isn’t responding press ‘control’ and the letter ‘c’ at the same time.

<p style="text-align: center;">
   <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise7.png"><img class="alignnone size-full wp-image-438" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise7.png" alt="Unix_exercise7" width="865" height="334" /></a>
</p>

**head**   
Display the first ten lines of a file

**tail**   
Display the last ten lines of a file  

Sometimes you may just want to view the text at the beginning or the end of a file, without having to display all of the file. The head and tail commands can be used to do this. You are now going to look at the beginning of S_typhi.embl.

`head S_typhi.embl [enter]`

To look at the end of S_typhi.embl type: `tail S_typhi.embl [enter]`
The number of lines that are displayed can be increased by adding extra arguments. To increase the number of lines viewed from 10 to 100 add the –100 argument to the command. For example to view the last 100 lines of S_typhi.embl type:

`tail -100 S_typhi.embl [enter]`

Do this for both head and tail commands. What type of information is at the beginning and end of the EMBL format file?

<p style="text-align: center;">
   <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise8.png"><img class="alignnone size-full wp-image-439" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise8.png" alt="Unix_exercise8" width="865" height="546" /></a>
</p>

**cat**  
Join files together. Having looked at the beginning and end of the S_typhi.embl file you should notice that in EMBL format files the annotation comes first, then the DNA sequence at the end. If you had two separate files containing the annotation and the DNA sequence, both in EMBL format, it is possible to concatenate or join the two together to make a single file like the S_typhi.embl file you have just looked at. The Unix command cat can be used to join two or more files into a single file. The order in which the files are joined is determined by the order in which they appear in the command line. For example, we have two separate files, MAL13P1.dna and MAL13P1.tab, that contain the DNA and annotation, respectively, from the *P. falciparum* genome. Return to the Module_Unix directory using the cd command:
`cd ../Module_Unix [enter]`

and type
`cat MAL13P1.tab MAL13P1.dna > MAL13P1.embl [enter]`  
MAL13P1.tab and MAL13P1.dna will be joined together and written to a file called MAL13P1.embl
The `>` symbol in the command line directs the output of the cat program to the designated file MAL13P1.embl

<p style="text-align: center;">
   <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise9.png"><img class="alignnone size-full wp-image-440" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise9.png" alt="Unix_exercise9" width="865" height="258" /></a>
</p>

**wc**  
Counts the lines, words or characters of files.
By typing the command line:
`ls | wc -l [enter]`

The above command uses wc to count the number of files that are listed by ls. The ‘-l’ option tells wc to return a count of the number of lines. The | symbol (known as the ‘pipe’ character) in the command line connects the two commands into a single operation for simplicity. You can connect as many commands as you want:

`ls | grep ".embl" | wc -l`
This command will list out all of the files in the current directory, then send the results to the grep command which searches for all filenames containing the ‘embl’, then sends the results to wc which counts the number of lines (which corresponds to the number of files).

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise10.png"><img class="alignnone size-full wp-image-441" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise10.png" alt="Unix_exercise10" width="865" height="212" /></a>
</p>

**grep**  
Searches a file for patterns. **grep** is a powerful tool to search for patterns in a file. In the examples below, we are going to use the file called Malaria.fasta that contains the set of *P. falciparum* chromosomes in FASTA format. A FASTA file has the following format:

> >Sequence Header
> CTAAACCTAAACCTAAACCCTGAACCCTAA...

Therefore if we want to get the sequence headers, we can extract the lines that match the ‘>’ symbol:

`grep ‘>’ Malaria.fasta [enter]`

By typing the command line:

`grep -B 1 -A 1 'aagtagggttca' Malaria.fasta [enter]`

This command will search for a nucliotide sequence and print 1 line before and after any match. It won’t find the pattern if it spans more than 1 line.

<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise11.png"><img class="alignnone size-full wp-image-442" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise11.png" alt="Unix_exercise11" width="865" height="536" /></a>
</p>

**find**
Finds files matching an expression. The find command is similar to ls but in many ways it is more powerful. It can be used to recursively search the directory tree for a specified path name, seeking files that match a given Boolean expression (a test which returns true or false)

`find . -name “*.embl”` This command will return the files which name has the .embl suffix.

`mkdir test_directory`
`find . -type d`

This command will return all the subdirectories contained in the current directory. These are just two basic examples but it is possible to search in many other ways: `-mtime` search files by modifying date `-atime` search files by last access date `-size` search files by file size `-user` search files by user they belong to.

#### Tips
You need to be careful with quoting when using wildcards!

The wildcard * symbol represents a string of any character and of any length.


<p style="text-align: center;">
  <a href="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise12.png"><img class="alignnone size-full wp-image-443" src="http://77.235.253.122/tutorials/wp-content/uploads/2015/06/Unix_exercise12.png" alt="Unix_exercise12" width="865" height="267" /></a>

For more information on Unix command see EMBNet UNIX Quick Guide.  

<h1 style="text-align: center;">
  End of the module
</h1>
# Introduction to Unix (continued)

In this part of the Unix tutorial, you will learn to download files, compress and decompress them, and combine commands.

## Download files

`wget` can be used to download files from the internet and store them.

`wget https://raw.githubusercontent.com/HadrienG/tutorials/master/LICENSE`

will download the file that is located at the above URL on the internet, and put it **in the current directory**. This is the license under which this course is released. Open it and read it if you like!

The `-O` option can be used to change the output file name.

`wget -O GNU_FDL.txt https://raw.githubusercontent.com/HadrienG/tutorials/master/LICENSE`

You can also use wget to download a file list using -i option and giving a text file containing file URLs. The following

```
cat > download-file-list.txt
url_1
url_2
url_3
url_4
[CTRL-C] (to exit cat)
```

`wget -i download-file-list.txt`

## Compressing and decompressing files

### Compressing files with gzip

gzip is a utility for compressing and decompressing individual files. To compress files, use:

`gzip filename`

The filename will be deleted and replaced by a compressed file called filename.gz. To reverse the compression process, use:

`gzip -d filename.gz`

Try it on the License you just downloaded!

### Tar archives

Quite often, you don't want to compress just one file, but rather a bunch of them, or a directory.

tar backs up entire directories and files as an archive. An archive is a file that contains other files plus information about them, such as their filename, owner, timestamps, and access permissions. tar does not perform any compression by default.

To create a gzipped disk file tar archive, use

`tar -czvf archivename filenames`

where archivename will usually have a .tar.gz extension

The c option means create, the v option means verbose (output filenames as they are archived), option f means file, and z means that the tar archive should be gzip compressed.

To list the contents of a gzipped tar archive, use

`tar -tzvf archivename`

To unpack files from a tar archive, use

`tar -xzvf archivename`

Try to archive the folder `Module_Unix` from the previous exercise!

You will notice a file called tutorials.tar.bz2 in your home directory. This is also a compressed archive, but compressed in the bzip format. Read the tar manual and find a way to decompress it.

Hint: you can read the manual for any command using `man`

`man tar`

### Redirection

Some commands give you an output to your screen, but you would have preferred it to go into another program or into a file. For those cases you have some redirection characters.

#### Output redirection

The output from a command normally intended for standard output (that is, your screen) can be easily diverted to a file instead. This capability is known as output redirection:

If the notation `> file` is appended to any command that normally writes its output to standard output, the output of that command will be written to file instead of your terminal.

I.e, the following who command:

`who > users.txt`

No output appears at the terminal. This is because the output has been redirected into the specified file.

`less users.txt`

Be careful, if a command has its output redirected to a file and the file already contains some data, that data will be lost. Consider this example:

`echo Hello > users.txt`

`less users.txt`

You can use the `>>` operator to append the output in an existing file as follows:

```
who > users.txt
echo "This goes at the end of the file" >> users.txt
```

`less users.txt`

#### Piping

You can connect two commands together so that the output from one program becomes the input of the next program. Two or more commands connected in this way form a pipe.

To make a pipe, put a vertical bar `|` on the command line between two commands.

Remember the command `grep`? We can pipe other commands to it, to refine searches per example:

`ls -l ngs_course_data | grep "Jan"`

will only give you the files and directories created in January.

Tip: There are various options you can use with the grep command, look at the manual!

Pipes are extremely useful to connect various bioinformatics software together. We'll use them extensively later.
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
# Introduction to UNIX (continued)

In the 4th and last module of your unix course, we'll how to write small programs, or scripts.

Shell scripts allow us to program commands in chains and have the system execute them as a scripted chain of events. They also allow for far more useful functions, such as command substitution. You can invoke a command, like date, and use it’s output as part of a file-naming scheme. You can automate backups and each copied file can have the current date appended to the end of its name. You can automate a bioinformatics analysis pipeline.

Before we begin our scripting tutorial, let’s cover some basic information. We’ll be using the bash shell, which most Linux distributions use natively. Bash is available for Mac OS users and Cygwin on Windows (which you are using with MobaXterm). Since it’s so universal, you should be able to script regardless of your platform.

At their core, scripts are just plain text files. You can use nano (or any other text editor) to write them.

## Permissions

Scripts are executed like programs. For this to happen, you need to have the proper permissions.
You can make the script executable for you by running the following command:

`chmod u+x my_script.sh`

by convention, bash script are saved with the .sh extension. Linux doesn't really care about file extension, but it is easier for the user to use the "proper" extensions!

## executing a script

You have to cd in the proper directory, then run the script like this:

`./my_script.sh`

To make things more convenient, you can place scripts in a “bin” folder in your home directory and add it to your path

`mkdir -p ~/bin`

More information on how to correctly modify your PATH [here](http://unix.stackexchange.com/a/26059)

## Getting started

As previously said, every script is a text file. Still, there are rules and conventions to follow in order of you file being recognized as a script

If you juste write a few command and try to execute it as is, with `./my_script`, it will not work. You can invoke `sh my_script`, but it is not very convenient.
`./` tries to find out which interpreter to use (e.g. which programming language and how to execute your script). It does so by looking at the first line:

The first line of your bash scripts should be:

`#!/bin/bash` or `#!/usr/bin/env bash`

The second version being better and more portable. Ask your teacher why!

This line will have the same syntax for every interpreted language. If you are programming in python:

`#!/usr/bin/env python`

### New line = new command

After the firstline, every line of your script will be a new command. Your first scripts will essentially be a succession of terminal commands. We'll learn about flow control (if, for, while, ...) later on.

### Comments

It is good practise to comment your scripts, i.e give some explanation of what is does, and explain a particularly arcane method that you wrote.

Comments start with a `#` and are snippets of texts that are ignored by the interpreter.

### Your first script

Let's start with a simple script, that copy files and append today's date to the end of the file name. We'll call it `datecp.sh`

In your `~/bin` folder:

```
touch datecp.sh
chmod u+x datecp.sh
```

and let's start writing our script

`nano datecp.sh`

```
#!/usr/bin/env bash

# this script will copy a file, appending the data and time to
# the end of the file name

```

Next, we need to declare a variable. A variable allows us to store and reuse information (characters, the date or the command `date`). Variables have a name, but can **expend** to their content when referenced if they contain a command.

Variables can hold strings and characers, like this:

`my_variable="hippopotamus"`

or a command. In bash, the correct way to store a command in a variable is within the syntax `$()`:

`variable=$(command –options arguments)`

Store the date and time in a variable. Test the date command first in your terminal, then when you got the right format, store it in a variable in your script.

It is generally bad practice to put spaces in file names in unix, so we'll want the following date format:

`date +%m_%d_%y-%H.%M.%S`

and for putting it into a variable:

date_formatted=$(date +%m_%d_%y-%H.%M.%S)

Your script now can print thedate without too much more coding:

```
#!/usr/bin/env bash

# this script will copy a file, appending the data and time to
# the end of the file name

date_formatted=$(date +%m_%d_%y-%H.%M.%S)
echo "This is the date and time:" $date_formatted

```

Now we need to add the copying part:

`cp –iv $1 $2.$date_formatted`

This will invoke the copy command, with two options: -i for asking for permission before overwriting a file, and -v for verbose.

You can also notice two variables: $1 and $2. When scripting in bash, a dollar sign ($) followed by a number will denote an argument of the script. For example in the following command:

`cp –iv a_file a_file_copy` the first argument ($1) is `a_file` and the second argument ($2) is `a_file_copy`

What our script will do is a simple copy of a file, but with adding the date to the end of the file name. Save it and try it out!

### Exercise

Write a script that backs itself up, that is, copies itself to a file named backup.sh.

Hint: Use the cat command and the appropriate positional parameter.



<!-- ## Conditionals

Conditionals let you decide whether to perform an action or not, this decision is taken by evaluating an expression.

### if / then / else -->

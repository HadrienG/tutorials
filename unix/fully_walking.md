# Introduction to Unix (continued)

In this part of the Unix tutorial, you will learn to download files, compress and decompress them, and combining commands

## Download files

`wget` can be used to download files from internet and store them. The following downloads and stores a file called to the current directory.

`wget https://raw.githubusercontent.com/HadrienG/tutorials/master/LICENSE`

will download the file that is located at the above URL on the internet, and put it **in the current directory**. This is the license under which this course is released. Open in and read it if you like!

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

The filename will be deleted and replaced by a compressed file called filename.gz To reverse the compression process, use:

`gzip -d filename.gz`

Try it on the License you just downloaded!

### Tar archives

Quite often, you don't want to compress just one file, but rather a bunch of them, or a directory.

tar backs up entire directories and files as an archive. An archive is a file that contains other files plus information about them, such as their filename, owner, timestamps, and access permissions. tar does not perform any compression by default.

To create a gzipped disk file tar archive, use

`tar -czvf archivename filenames`

where archivename will usually have a .tar .gz extension

The c option means create, the v option means verbose (output filenames as they are archived), and option f means file.

To list the contents of a gzipped tar archive, use

`tar -tzvf archivename`

To unpack files from a tar archive, use

`tar -xzvf archivename`

Try to archive the folder `Module_Unix` from the previous exercise!

You will notice a file called tutorials.tar.bz2 in your home directory. This is also a compressed archive, but compressed in the bzip format. Read the tar manual and find a way to decompress it

Hint: you can read the manual for any command using `man`

`man tar`

### Redirection

Some commands give you an output to your screen, but you would have preferred it to go into another program or into a file. For those cases you have some redirection characters.

#### Output redirection

The output from a command normally intended for standard output (that is, your screen) can be easily diverted to a file instead. This capability is known as output redirection:

If the notation `> file` is appended to any command that normally writes its output to standard output, the output of that command will be written to file instead of your terminal

I.e, the following who command:

`who > users.txt`

No output appears at the terminal. This is because the output has been redirected  into the specified file.

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

will only give you the files and directories created in January

Tip: There are various options you can use with the grep command, look at the manual!

Pipes are extremely useful to connect various bioinformatics softwares together. We'll use them extensively later.

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

$ tar -xzvf archivename

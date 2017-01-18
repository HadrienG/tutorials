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

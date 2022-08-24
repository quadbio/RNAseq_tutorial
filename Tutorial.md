# Tutorial for bulk RNA-seq data preprocessing and analysis
#### Compiled by Zhisong He
#### Updated on 24 Aug 2022
### Table of Content
  * [Introduction](#introduction)
  * [Preparation](#preparation)
    * [1-1. Linux and Bash, your important buddies for RNA-seq data analysis](#1-1-linux-and-bash-your-important-buddies-for-rna-seq-data-analysis)
  * Preprocessing of RNA-seq data
  * Analyze and compare RNA-seq data

## Introduction
Nowadays in biology, the RNA sequencing technology becomes so commonly used that it may not even need much of the introduction. Based on the high-throughput DNA sequencing techniques, the RNA sequencing technology, or RNA-seq in short, measures the presence and quantity of RNA in a biological sample by firstly converting the captured RNA transcripts to cDNA via reverse transcription. Unlike the other techniques which also quantify the transcript abundance (e.g. qPCR), RNA-seq can unbiasedly quantify RNA abundance in the transcriptome level at the same time, and in addition allows the discovery of novel transcripts in a sample.

On the other hand, no RNA-seq experiment would automatically inform any of those biological insights directly. What an RNA-seq experiment gives you is a bunch of "reads", which are the short fragmented sequences that are captured, fragmented, reverse transcribed and measured. One would have to process the data to quantify the abundance of transcripts from every gene, and then further analyze the resulted data. This is what this tutorial is supposed to tell you.

In the following sections, we will go through the RNA-seq data preprocessing including quality control, mapping and quantification, and then the data analysis such as sample/condition comparison, identification of differentially expressed genes among conditions, and the possible following analysis to interpret the identified differences. Before that, I will also mention some more basic concepts and operations such as how to use the command line interface in a Unix-like system like Linux, how to install the software needed for the procedure, and how to download public RNA-seq data from the SRA data repository.

## Preparation
In this section, I will mention the stuffs before we even start preprocessing the RNA-seq data, so that you can make sure that you and your computers are both ready for the way afterwards. These stuffs include the followings:
1. A brief introduction and simple HOW-TO of [Bash (the Bourne Again SHell)](https://www.gnu.org/software/bash/), the most commonly used command line interface in Linux
2. A brief introduction of conda and how to use it to install software needed for the data preprocessing
3. A brief introduction of SRA and how to retrieve the public RNA-seq data from there

### 1-1. Linux and Bash, your important buddies for RNA-seq data analysis
I'm pretty sure you know computer and use computer (otherwise you won't be able to see this). However, you may not ever use Linux OS (Linux for short) even though you may have heard of the name somewhere. Meanwhile, you must be very familiar with Microsoft Windows and/or Apple MacOS. These three falls into the same concept: the operating system (OS) of a computer. An OS is probably the most important software on a computer, which manages computer hardware, software resource and provides services for other programs. It runs at the lowest level, and everything else you use on the computer relies on it. On desktop computers or laptops, Microsoft Windows and Apple MacOS are the two most commonly used OS. However, this is not the case for other computers such as computing servers and high-performance clusters (HPC), which are important foundations of the nowaday data science, and of course, bioinformatics. For those computers, Linux is of much higher preference, for its computing efficiency, flexibility, stability, and security. Indeed, more than 90% of the world's fastest supercomputers run on Linux, and because of that, people relying on high-performance computing develop lots of tools that can be quite easily set up in Linux but not necessarily in the other two OS, which also contributes to the bias of choice.

As one of the fields that require high-performancing and large-resource computing, bioinformatics and computational biology also heavily uses Linux. Actually, <ins>**many software needed for the RNA-seq data preprocessing are only available in Linux**</ins> (and other Unix-like systems, but won't be further mentioned in this tutorial). Therefore, some basic understanding of how to use Linux is critical.

<span style="font-size:0.8em">*P.S. Unix is a family of multitasking, multiuser computer operating systems that derive from the original AT&T Unix, and are characterized by a modular design (Unix philosophy). Different components work together, with the kernel being the center to provide the most basic management and support. Unix is not free, but it inspired the development of many free Unix and Unix-like OS, and importantly, the GNU (meaning GNU's Not Unix) project and the Linux kernel which later on further derived into different Linux distributions (like different branches or versions of Linux). Among them, Red Hat Enterprise Linux (RHEL), Fedora, CentOS, SUSE Linux Enterprise, openSUSE, Debian, Ubuntu, Linux Mint and Arch Linux are among the most famous ones.*</span>

<span style="font-size:0.8em">*P.S.2. Actually the Apple MacOS (since version 10.5 with 10.7 Lion as the exception) is a UNIX 03-compliant OS certified by The Open Group, so although not a Unix-like (meaning look like Unix but is not been certified, e.g. Linux) but a UNIX. This is also the reason why MacOS in prefered in relative to Windows if you have to use your desktop/laptop for the task.*</span>

From our side as the computer end users, we don't need to care too much about how different OS work at the low level. However, we need to interact with the OS, for instance, to ask it to open a software or whatever, and that varies a lot from one OS to another. Windows and MacOS look very different and we know that. When switching from one to the other for the first time, most of us would need quite a long time to get used to everything. However, that usually wouldn't be too difficult as both OS have pretty and straightforward graphical design for you to interact with the system. This is called a graphic user interface (GUI). Linux also have GUI, and different Linux distributions have different ones. However, what is different from Windows and MacOS is that the GUI is not essential for Linux. For many computers running on Linux, especially those for high-performance computing, the GUI components are not even installed. Instead, people using those computers rely on CLI, which is short for command-line interface.

<p align="center">
<img src="img/GUI_vs_CLI.jpg"/>
<br/><span style="font-size:0.8em"><i>Image adapted from https://www.onetekno.my.id/2021/12/perbedaan-gui-dan-cli.html.</i></span>
</p>

CLI has a much steeper learning curve than GUI, and that's exactly the reason why GUI was developed in response to the complain to CLI. So why are we still using CLI heavily?That was also a question I had at the beginning, but after using CLI for a while and finally getting into it, I realized at least several reasons.

1. <ins>**The fancy graphic view comes with cost**</ins>. Computing and displaying all the visual changes on the screen needs resources (the computing unit CPU and/or GPU, memory, and also the storage). When doing high-performance computing which requires a lot of resource, however, we would for sure want to maximize the resource we can use for the real computation we want to do rather than just to display a mouse cursor moving around. Don't underestimate how much this would have needed, especially keep in mind that many computing servers and HPCs are not just been used by just one person at a time, and one person may have many things running at the same time. With a GUI for every user at least would use tremendous amount of computer resource.
2. While the GUI can only do what the developers implemented explicitly, <ins>**the CLI is easy to program to work out something more complicated or repetitive**</ins>. For instance, renaming 1000 files with a similar manner (e.g. to rename everything from A\*\*\*\* to B\*\*\*\* while keeping the \*\*\*\* part the same) can be quite easily done in CLI using the basic rename command together with basic programming skills in several lines of commands, but you would probably need to find a tool specifically implemented for batch renaming to do that in GUI. One can also combine different steps which use different tools together quite easily in CLI, much easier than using GUI. Such kinds of operations are very common when dealing with large data (e.g. preprocessing a RNA-seq data set with 10 samples using the same pipeline).
3. <ins>**Developing GUI for software needs a lot of effort**</ins>. Designing and implementing a nice GUI for a program can take as much time as, if not more than, implementing the functional part of the program. For tools which are mostly used by CLI-familiar users, it doesn't make sense from the developer's perspective to implement the GUI. If all the tools running on the OS only use CLI, the GUI of the OS is not really that useful anymore.

Learning and getting used to CLI would be pretty tough at the beginning, but it is far from doing something impossible, and you would probably like it once you get familiar with the commands and learn a bit of the simple programming skills. It takes time and effort of course, but stay patient and believe in yourself that you can do it! 

Before going to the HOW-TO part, do keep in mind that there are different CLIs for different OS (yes, Windows and MacOS both have CLI as well). And there are different CLIs implemented for even the same OS. For Windows, there are the Command Prompt which emulates many of the command line abilities available in MS-DOS (Microsoft Disk Operating System, the old OS with only CLI by Microsoft), and the new PowerShell which provides extended capacities of the Command Prompt. For UNIX and Unix-like OS including Linux and MacOS, there are different types of Shell including C shell (csh), Bourne Again shell (Bash), and Z shell (zsh), an extension of Bash. In most of the time, Bash is the default CLI for a Linux regardless different distributions, while zsh is currently the default shell for MacOS.

Here let's focus on Bash, and go through a few commonly used commands in Bash. Please see it as a start of getting into using Bash in Linux, not the end. If you would like to dedicate yourself a bit into data analysis, knowing more than the most basic commands would be very useful and important. For instance, the scripting function is extremely useful to wrap up different operations into one pipeline and apply to multiple objects (like files), but this won't be covered here (otherwise this would become a Bash tutorial than RNA-seq data analysis). There are many great books introducing Bash commands and scripting, so as many resources available online. You can quite easily get the information online.

<span style="font-size:0.8em">*P.S. The term "shell" here means a computer program that presents a CLI that allows you to control your computer using commands entered with a keyboard instead of controlling GUIs with a mouse/keyboard/touchscreen combination.*</span>

<span style="font-size:0.8em">*P.S.2 The term prompt refers to what you see when the CLI is waiting for your command. In Bash it is`$` by default but customizable.*</span>

<span style="font-size:0.8em">*P.S.3 The Windows CLIs has a lot of differences compared to the shells used in Linux, and we won't talk about it here.*</span>

Here are some of the most commonly used Bash commands, most of which are pre-installed in any commonly used Linux distribution:

| Command | Abbrev. | Function |
|---------|---------|----------|
|`ls`|list|list files, by default files in the current working folder, but could also be those at a given location or those with their paths following a certain pattern.|
|`cd`|change directory|change the working folder|
|`pwd`|print working directory|print the current working directory|
|`touch`|touch|update the access and modification times of the given file to the current time. Usually it is used to create an empty file by providing a file name that points to no file at the moment|
|`cp`|copy|copy a file or files to a different path (could be a file named the same or differently in another folder, or a different file name at the same folder)|
|`rm`|remove|delete the given file(s)|
|`mv`|move|move a file or files to a different path, or to rename it/them. Similar to firstly do `cp` and then `rm`|
|`mkdir`|make directory|create a new and empty folder at the current working folder or the given path|
|`rmdir`|remove directory|remove the given folder(s). Note that this only applies to empty folders|
|`chmod`|change mode|manage the access permissions of a file|
|`cat`|concatenate|output the content of one or multiple files. When multiple files are given, the contents are concatenated (that's what the name is for)|
|`less`| |view the contents of a text file, one screen at a time. One can use the PgUp and PgDn buttom to go forward or backward|
|`echo`|echo|print the given strings|
|`gzip`|GNU zip|compress or decompress a single file (directory not allowed) based on the gzip file format|
|`tar`|tape archive|collect many files into one archived file (tarball), or to release the files from the given tarball|
|`which`|which|print the location of the given executable|
|`wget`|www get|download files via HTTP/HTTPS/FTP from Internet, given the link to the files|
|`ssh`|secure shell|remotely login and access another computer using the SSH network communication protocol|
|`scp`|secure copy|transfer files between the local computer and the remote |
|`wc`|word count|count the number of characters, words and lines in a file or files|
|`top`|table of processes|display information about CPU and memory utilzation in total and per process|
|`screen`|screen|open one or more virtual screens which can be stayed in the background while keeping commands going on even if the virtual screen is disconnected|
|`tmux`|terminal multiplexer|very similar to `screen`, just another option|
|`nano`||an easy-to-use text editor|
|`vim`|vi improved|a commonly used modal text editor, which operates in either insert mode (where typed text becomes part of the document) or command mode (where keystrokes are interpreted as commands that control the edit session)|

<span style="font-size:0.8em">*P.S. Any command you can run in Bash is essentially a program. On the other hand, any program which is executable in Bash, including those implemented separately from Bash, can be also seen as a command there.*</span>

All those commands have different options and arguments. In most of the time, they are arranged in a way like following:
```
$ <command> [options] [arguments]
```

Here, `arguments` represent one or multiple things (e.g. the path to a file, the link to a file on Internet) as the given input to the command, and it is often required unless the command provides a default value (e.g. `ls` has the default argument `.` which means the current working directory). Meanwhile, `options` are the parameters specific to commands which change the behavior of the command. For instance, the `ls` command has many `options`, such as `-a` for displaying all files including the hidden ones, and `-l` for showing the detailed information of all files with each line per file. Note that it is possible that certain `options` changes the command behavior so that no `argument` is anymore expected. The common example is the `-h` or `--help` option which usually asks the command to show the brief or detailed description of possible `options` and the expected `arguments`. For instance, running `ls --help` shows the following (it is too long so only the beginning is shown here):

```
$ ls --help
Usage: ls [OPTION]... [FILE]...
List information about the FILEs (the current directory by default).
Sort entries alphabetically if none of -cftuvSUX nor --sort is specified.

Mandatory arguments to long options are mandatory for short options too.
  -a, --all                  do not ignore entries starting with .
  -A, --almost-all           do not list implied . and ..
      --author               with -l, print the author of each file
  -b, --escape               print C-style escapes for nongraphic characters
      --block-size=SIZE      scale sizes by SIZE before printing them; e.g.,
                               '--block-size=M' prints sizes in units of
                               1,048,576 bytes; see SIZE format below
  -B, --ignore-backups       do not list implied entries ending with ~
  -c                         with -lt: sort by, and show, ctime (time of last
                               modification of file status information);
                               with -l: show ctime and sort by name;
                               otherwise: sort by ctime, newest first
  -C                         list entries by columns
      --color[=WHEN]         colorize the output; WHEN can be 'never', 'auto',
                               or 'always' (the default); more info below
  -d, --directory            list directories themselves, not their contents
```

Another way to view the detailed description of a command is to use the `man` command, which displays the manual of the given command if the manual is available. For instance, running `man ls` shows the system manual of the `ls` command. And of course, you can run `man man` to see the manual of the `man` command.

<span style="font-size:0.8em">*P.S. In the manual page by `man`, use PgUp and PgDn to scroll, and q to quit. The same way also applies when you use the `less` command to view a text file.*</span>

Some simple examples of full command lines in Bash:

|Command|Function|
|--------|--------|
|`cd folder`|change the working directory to the `folder` directory at the current working directory|
|`cd /tmp`| change the working directory to `/tmp`|
|`ls`|display files in the current working directory|
|`ls -l`|display files in the current working directory with detailed information|
|`ls -al`|display all files, including the hidden ones, in the current working directory with detailed information|
|`ls -hl`|display files in the current working directory, with the file size displayed in a human-readable format|
|`ls /tmp`|display files in the `/tmp` directory|
|`cp file1 file2`|make a copy of the file `file1` in the current working directory as the `file2` file in the same folder|
|`cp file1 /tmp`|make a copy of the file `file1` in the current working directory in the `/tmp` folder with the same name|
|`cp file1 /tmp/file2`|make a copy of the file `file1` in the current working directory as the `file2` file the `/tmp` folder|
|`cp -r dir1 /tmp/dir2`|make a recursive copy of the folder `dir1` as the folder `dir2` in the `/tmp` folder|
|`rm file1`|delete the file `file1` in the current working directory|
|`rm -r dir1`|delete the folder `dir1` and all its contents in the current working directory|
|`gzip file1`|compress `file1` to `file1.gz`|
|`gzip -d file1.gz`|decompress `file1.gz` to `file1`|

Before ending this section, there are several extra things that I think should worth a mentioning.

One thing is the file structure in Linux and other Unix-like and UNIX OS (such as MacOS). Files are always organized in a hierarchical manner, with the top level being the root `/`. Under `/` there are different directories, each is specifically for one purpose. For instance, `/home` contains all the user home folders, `/usr` contains most of the commands and software accessible to all users, `/tmp` is the default temporary folder to store temporary files. To represent a file under a certain folder, concatenate the file name with its higher hierarchy by `/`, e.g. `/home/user1`, and this is called the absolute path to a file. There is also the relative path, meaning the path of a file in relative to another file (in most of the time the current working directory). This would need to include `.` for the current working directory, and `..` to represent the previous hierarchy (or the parent directory) of the current working directory. For instance, assuming the current working directory is `/home/user1`, then `../user2` means the file `/home/user2`. By the way, the terms "directory" and "folder" are somehow equivalent here.

You may have also noticed that I use the term "file" to represent not only the literally files, but also directories. Indeed in Linux, everything is represented or can be seen as a file. This not only include directories, but also hardware devices being attached and recognized. They are all seen and managed as special files with some special behaviors.

[pipe]

[scripting]



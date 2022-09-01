# Tutorial for bulk RNA-seq data preprocessing and analysis
#### Compiled by Zhisong He
#### Updated on 29 Aug 2022
### Table of Content
  * [Introduction](#introduction)
  * [Preparation](#preparation)
    * [1-1. Linux and Bash, your important buddies for RNA-seq data analysis](#1-1-linux-and-bash-your-important-buddies-for-rna-seq-data-analysis)
    * [1-2. Access the computing server](#1-2-access-the-computing-server)
    * [1-3. Install the required tools with the help from conda](#1-3-install-the-required-tools-with-the-help-from-conda)
    * [1-4. Get the public RNA-seq data from SRA](#1-4-get-the-public-rna-seq-data-from-sra)
  * [Preprocessing of RNA-seq data](#preprocessing-of-rna-seq-data)
    * [2-1 Quality control of RNA-seq data](#2-1-quality-control-of-rna-seq-data)
    * [2-2 Read mapping/pseudomapping and quantification](#2-2-read-mappingpseudomapping-and-quantification)
      * [2-2-1 Read mapping via STAR and data quantification](#2-2-1-read-mapping-via-star-and-data-quantification)
  * Analyze and compare RNA-seq data

## Introduction
Nowadays in biology, the RNA sequencing technology becomes so commonly used that it may not even need much of the introduction. Based on the high-throughput DNA sequencing techniques, the RNA sequencing technology, or RNA-seq in short, measures the presence and quantity of RNA in a biological sample by firstly converting the captured RNA transcripts to cDNA via reverse transcription. Unlike the other techniques which also quantify the transcript abundance (e.g. qPCR), RNA-seq can unbiasedly quantify RNA abundance in the transcriptome level at the same time, and in addition allows the discovery of novel transcripts in a sample.

On the other hand, no RNA-seq experiment would automatically inform any of those biological insights directly. What an RNA-seq experiment gives you is a bunch of "reads", which are the short fragmented sequences that are captured, fragmented, reverse transcribed and measured. One would have to process the data to quantify the abundance of transcripts from every gene, and then further analyze the resulted data. This is what this tutorial is supposed to tell you.

In the following sections, we will go through the RNA-seq data preprocessing including quality control, mapping and quantification, and then the data analysis such as sample/condition comparison, identification of differentially expressed genes among conditions, and the possible following analysis to interpret the identified differences. Before that, I will also mention some more basic concepts and operations such as how to use the command line interface in a Unix-like system like Linux, how to install the software needed for the procedure, and how to download public RNA-seq data from the SRA data repository.

## Preparation
In this section, I will mention the stuffs before we even start preprocessing the RNA-seq data, so that you can make sure that you and your computers are both ready for the way afterwards. These stuffs include the followings:
1. Linux and [Bash (the Bourne Again SHell)](https://www.gnu.org/software/bash/), the most commonly used command line interface in Linux
2. How to access a computing server via SSH
3. Conda and how to use it to install software needed for the data preprocessing
4. SRA and how to retrieve the public RNA-seq data from there

### 1-1. Linux and Bash, your important buddies for RNA-seq data analysis
I'm pretty sure you know computer and use computer (otherwise you won't be able to see this). However, you may not ever use Linux OS (Linux for short) even though you may have heard of the name somewhere. Meanwhile, you must be very familiar with Microsoft Windows and/or Apple macOS. These three falls into the same concept: the operating system (OS) of a computer. An OS is probably the most important software on a computer, which manages computer hardware, software resource and provides services for other programs. It runs at the lowest level, and everything else you use on the computer relies on it. On desktop computers or laptops, Microsoft Windows and Apple macOS are the two most commonly used OS. However, this is not the case for other computers such as computing servers and high-performance clusters (HPC), which are important foundations of the nowaday data science, and of course, bioinformatics. For those computers, Linux is of much higher preference, for its computing efficiency, flexibility, stability, and security. Indeed, more than 90% of the world's fastest supercomputers run on Linux, and because of that, people relying on high-performance computing develop lots of tools that can be quite easily set up in Linux but not necessarily in the other two OS, which also contributes to the bias of choice.

As one of the fields that require high-performancing and large-resource computing, bioinformatics and computational biology also heavily uses Linux. Actually, <ins>**many software needed for the RNA-seq data preprocessing are only available in Linux**</ins> (and other Unix-like systems, but won't be further mentioned in this tutorial). Therefore, some basic understanding of how to use Linux is critical.

>**NOTE**
>* Unix is a family of multitasking, multiuser computer operating systems that derive from the original AT&T Unix, and are characterized by a modular design (Unix philosophy). Different components work together, with the kernel being the center to provide the most basic management and support. Unix is not free, but it inspired the development of many free Unix and Unix-like OS, and importantly, the GNU (meaning GNU's Not Unix) project and the Linux kernel which later on further derived into different Linux distributions (like different branches or versions of Linux). Among them, Red Hat Enterprise Linux (RHEL), Fedora, CentOS, SUSE Linux Enterprise, openSUSE, Debian, Ubuntu, Linux Mint and Arch Linux are among the most famous ones.
>* Actually the Apple macOS (since version 10.5 with 10.7 Lion as the exception) is a UNIX 03-compliant OS certified by The Open Group, so although not a Unix-like (meaning look like Unix but is not been certified, e.g. Linux) but a UNIX. This is also the reason why macOS in prefered in relative to Windows if you have to use your desktop/laptop for the task.

From our side as the computer end users, we don't need to care too much about how different OS work at the low level. However, we need to interact with the OS, for instance, to ask it to open a software or whatever, and that varies a lot from one OS to another. Windows and macOS look very different and we know that. When switching from one to the other for the first time, most of us would need quite a long time to get used to everything. However, that usually wouldn't be too difficult as both OS have pretty and straightforward graphical design for you to interact with the system. This is called a graphic user interface (GUI). Linux also have GUI, and different Linux distributions have different ones. However, what is different from Windows and macOS is that the GUI is not essential for Linux. For many computers running on Linux, especially those for high-performance computing, the GUI components are not even installed. Instead, people using those computers rely on CLI, which is short for command-line interface.

<p align="center">
<img src="img/GUI_vs_CLI.jpg"/>
<br/><sub><i>Image adapted from https://www.onetekno.my.id/2021/12/perbedaan-gui-dan-cli.html.</i></sub>
</p>

CLI has a much steeper learning curve than GUI, and that's exactly the reason why GUI was developed in response to the complain to CLI. So why are we still using CLI heavily?That was also a question I had at the beginning, but after using CLI for a while and finally getting into it, I realized at least several reasons.

1. <ins>**The fancy graphic view comes with cost**</ins>. Computing and displaying all the visual changes on the screen needs resources (the computing unit CPU and/or GPU, memory, and also the storage). When doing high-performance computing which requires a lot of resource, however, we would for sure want to maximize the resource we can use for the real computation we want to do rather than just to display a mouse cursor moving around. Don't underestimate how much this would have needed, especially keep in mind that many computing servers and HPCs are not just been used by just one person at a time, and one person may have many things running at the same time. With a GUI for every user at least would use tremendous amount of computer resource.
2. While the GUI can only do what the developers implemented explicitly, <ins>**the CLI is easy to program to work out something more complicated or repetitive**</ins>. For instance, renaming 1000 files with a similar manner (e.g. to rename everything from A\*\*\*\* to B\*\*\*\* while keeping the \*\*\*\* part the same) can be quite easily done in CLI using the basic rename command together with basic programming skills in several lines of commands, but you would probably need to find a tool specifically implemented for batch renaming to do that in GUI. One can also combine different steps which use different tools together quite easily in CLI, much easier than using GUI. Such kinds of operations are very common when dealing with large data (e.g. preprocessing a RNA-seq data set with 10 samples using the same pipeline).
3. <ins>**Developing GUI for software needs a lot of effort**</ins>. Designing and implementing a nice GUI for a program can take as much time as, if not more than, implementing the functional part of the program. For tools which are mostly used by CLI-familiar users, it doesn't make sense from the developer's perspective to implement the GUI. If all the tools running on the OS only use CLI, the GUI of the OS is not really that useful anymore.

Learning and getting used to CLI would be pretty tough at the beginning, but it is far from doing something impossible, and you would probably like it once you get familiar with the commands and learn a bit of the simple programming skills. It takes time and effort of course, but stay patient and believe in yourself that you can do it! 

Before going to the HOW-TO part, do keep in mind that there are different CLIs for different OS (yes, Windows and macOS both have CLI as well). And there are different CLIs implemented for even the same OS. For Windows, there are the Command Prompt which emulates many of the command line abilities available in MS-DOS (Microsoft Disk Operating System, the old OS with only CLI by Microsoft), and the new PowerShell which provides extended capacities of the Command Prompt. For UNIX and Unix-like OS including Linux and macOS, there are different types of Shell including C shell (csh), Bourne Again shell (Bash), and Z shell (zsh), an extension of Bash. In most of the time, Bash is the default CLI for a Linux regardless different distributions, while zsh is currently the default shell for macOS.

Here let's focus on Bash, and go through a few commonly used commands in Bash. Please see it as a start of getting into using Bash in Linux, not the end. If you would like to dedicate yourself a bit into data analysis, knowing more than the most basic commands would be very useful and important. For instance, the scripting function is extremely useful to wrap up different operations into one pipeline and apply to multiple objects (like files), but this won't be covered here (otherwise this would become a Bash tutorial than RNA-seq data analysis). There are many great books introducing Bash commands and scripting, so as many resources available online. You can quite easily get the information online.

>**NOTE**
>* The term "shell" here means a computer program that presents a CLI that allows you to control your computer using commands entered with a keyboard instead of controlling GUIs with a mouse/keyboard/touchscreen combination.
>* The term prompt refers to what you see when the CLI is waiting for your command. In Bash it is`$` by default but customizable.
>* The Windows CLIs has a lot of differences compared to the shells used in Linux, and we won't talk about it here.

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

>**NOTE**
>Any command you can run in Bash is essentially a program. On the other hand, any program which is executable in Bash, including those implemented separately from Bash, can be also seen as a command there.

All those commands have different options and arguments. In most of the time, they are arranged in a way like following:
```console
$ <command> [options] [arguments]
```

Here, `arguments` represent one or multiple things (e.g. the path to a file, the link to a file on Internet) as the given input to the command, and it is often required unless the command provides a default value (e.g. `ls` has the default argument `.` which means the current working directory). Meanwhile, `options` are the parameters specific to commands which change the behavior of the command. For instance, the `ls` command has many `options`, such as `-a` for displaying all files including the hidden ones, and `-l` for showing the detailed information of all files with each line per file. Note that it is possible that certain `options` changes the command behavior so that no `argument` is anymore expected. The common example is the `-h` or `--help` option which usually asks the command to show the brief or detailed description of possible `options` and the expected `arguments`. For instance, running `ls --help` shows the following (it is too long so only the beginning is shown here):

```console
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

>**NOTE**
>In the manual page by `man`, use PgUp and PgDn to scroll, and q to quit. The same way also applies when you use the `less` command to view a text file.

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

#### About the file system in Linux
One thing is the file structure in Linux and other Unix-like and UNIX OS (such as macOS). Files are always organized in a hierarchical manner, with the top level being the root `/`. Under `/` there are different directories, each is specifically for one purpose. For instance, `/home` contains all the user home folders, `/usr` contains most of the commands and software accessible to all users, `/tmp` is the default temporary folder to store temporary files. To represent a file under a certain folder, concatenate the file name with its higher hierarchy by `/`, e.g. `/home/user1`, and this is called the absolute path to a file. There is also the relative path, meaning the path of a file in relative to another file (in most of the time the current working directory). This would need to include `.` for the current working directory, and `..` to represent the previous hierarchy (or the parent directory) of the current working directory. For instance, assuming the current working directory is `/home/user1`, then `../user2` means the file `/home/user2`. By the way, the terms "directory" and "folder" are somehow equivalent here.

You may have also noticed that I use the term "file" to represent not only the literally files, but also directories. Indeed in Linux, everything is represented or can be seen as a file. This not only include directories, but also hardware devices being attached and recognized. They are all seen and managed as special files with some special behaviors.

#### About the filename patterns
And when talking about the file names, one doesn't always need to provide the exact file name. Alternatively, one can represent a file, or a set of files, using the so-call glob patterns that specify a set of filenames following similar patterns with wildcard characters. Those wildcard characters are like placeholders with just one character. The commonly used ones include
 * the asterisk character `*` (or "star"), which represents zero or more of any characters
 * the question mark `?`, which represents exactly one of any character
 * the enclosed square branckets `[...]` with a set of characters in between (e.g. `[ab]`), which represent a single character that matches within the set
 * the exclamation mark within the enclosed square branckets as the first character `[!...]` (e.g. `[!ab]`), which represent a single character that matches with any character that's not in the set

Using those wildcard characters we can easily represent many different filenames that share a similar pattern. For instance, `a*` means any file starting with "a", `*a*` means any file with an "a", `*a` means any file ending with "a", `a*a` means any file starting with "a" and ending with "a", `a?a` means any filename with three characters in which the first and the last one being "a", `a[abc]a` means the set of `aaa`, `aba` and `aca`. And just to clarify, different wildcard characters can be used together, so we can have something like this `a[abcd]??[0123]*.txt`. The glob pattern can be used as `arguments` of many commands that expect one or multiple filenames, e.g. `ls` (`ls -l *.csv`, list all files ends with ".csv") and `cp` (`cp *.txt /tmp/`, copy all files end with ".csv" to the `/tmp` directory)

#### About redirection and piping
After running a certain command in Bash, you usually see the output printed on the screen. However, sometimes you would want to store all those results to a text file. This is possible by using  the redirection function in Bash, which is activated by using `>` after a command:

```console
<command> [options] [arguments] > [file]
```

For instance, you can save the output of `ls -l` to the file "ls.txt" by doing `ls -l > ls.txt`. This doesn't seem to be very useful, but the redirection function could become very useful in some scenarios. For instance, the `gzip -d` function decompresses a .gz file and afterwards the .gz file will be automatically removed. The `gzip` command has another option (`-c`) to print the compression/decompression results to the screen. So if you want to decompress a .gz file and save the decompressed context to a new file while keeping the .gz file, you can use the redirection function as

```console
gzip -cd example.txt.gz > example.txt
```

>**NOTE**
>There are two types of output to the screen, one is called "stdout" (standard output) while the other one called "stderr" (standard error). There is also the "stdin" (standard input) which is usually the input from the keyboard. Formally they are called standard I/O streams, I/O means input/output, and streams here represent the flow of information. For the two standard output streams, stdout is usually the real output while stderr is usually for verbose, warning, or error message. The simple `>` only save stdout to the file. If you want to save the stderr, `2>` instead of `>` should be used, where "2" represents stderr. And by the way, as you may have guessed, "1" represent stdout, so you can also use `1>` which is actually the same as using `>` directly.

Besides, there is another important feature of Bash, the pipes, indicating by `|`. It is a bit similar to redirection, but instead of saving the output to a file, it directly use the output of one command (which would be printed to the screen if you just run that command directly) as the input of the next command. This could be extremely useful to combine multiple commands for some complicated operations without the need to generate any intermediates.

```console
<command1> | <command2>
```

For example, the `ls` command has the option `-1` to print one file per line. Meanwhile, the `wc` command has the option `-l` to only output the number of lines in the given file or the piping input. We can therefore combine them two using the pipes

```console
ls -1 | wc -l
```

In this combination, the output of `ls -1` becomes the input of `wc -l`, so the number of lines in the `ls -1` is printed. As the line number of `ls -1` is the same as the number of files in the current working directory, the final output actually tells you how many files there are in the current folder.

>**NOTE**
>One can do multiple piping to build a pipeline. Of course, that would require that every command being used in the pipeline supports the use of stdin as the input and can output their results to stdout so that it can be piped into the next command.


#### About Bash scripting
Do keep in mind, that Bash can do much more complicated things than what have been mentioned above. It supports scripting, which allows a series commands being put together, plus the additional logic operators, for loops, conditional statements and so on. So you can actually see Bash as a programming language. Indeed, this is probably the most critical and valuable part of Bash. If you want to be an expert on Bash, this is what you have to learn.

Specifically for the main topic of this tutorial, to preprocess and analyze RNA-seq data, this would also be very helpful even with some simple knowledge on it. For example, assuming you have the data of 20 different RNA-seq samples and you know the preprocessing pipeline, you can of course apply the pipeline to each sample one by one, to manually start the next one after seeing the previous one being done, but this is clearly not the optimal way as you don't want to look at the screen 24\*7. And no need to mention when you have 200 samples instead of 20. This can be however easily managed with looping through the samples, so that every time one sample is being processed and the next one would be automatically started when the previous one is finished.

However, there are so many stuffs one would have to talk about the scripting in Bash, as many as introducing a programming language, it won't be covered further here. Meanwhile, there are lots of great materials one can find online (for instance, this [cheatsheet](https://devhints.io/bash)), as well as books. If possible, I would really encourage you to check and get a bit into it, and very likely you won't regret.

### 1-2. Access the computing server 
It is most likely that you are using your own computer right now looking at this tutorial, and you might be now thinking to use it for the RNA-seq data preprocessing and analysis, at least to go through the remaining part of the tutorial. This may not be a good idea, as your own computer, which is very much likely to be an ordinary desktop or laptop. However, while it could have been powerful enough for daily usage for sure and maybe also to play fancy AAA games without frame drops, it may not be enough for going through every step in the tutorial. Different tasks need different resources. A powerful graphics card is critical for gaming and on top of that, an Intel Core-i9 12th gen with 16 cores plus 16GB RAM (memory) would have made your computer powerful enough for most games. On the other hand, for most of the preprocessing and analysis of RNA-seq data, at least those will be mentioned in this tutorial, no graphic card would be needed at all, while RAM would become a main issue, that 64GB or more would be expected to get everything done; 16 cores CPU would be probably enough in most of the time, but sometimes during the analysis you may want to use large-scale parallelization to speed things up, and start to hope for availability of 100 cores.

Therefore, we rarely use personal computers for the data analysis. Instead, we use high performance computing servers or clusters. They are machines with large number of CPU cores and memories, and usually shared by many users who can all access and use the machine at the same time. And of course, this "access" would no longer be the physical access sitting in front of its screen and type directly on its keyboard. Actually many of those machine has no screen and keyboard (and of course no mouse as well) attached. Then how shall we use them?

The solution has actually mentioned somewhere above. Basically all of those machines supports remote access via the SSH (Secure Shell) protocol. Under that protocol, the server provides the service as a server, and we the users are the clients that try to request the services. The way to do this request varies depending on the OS you use at your personal computer.

#### If you are a macOS/Linux/other UNIX/Unix-like OS user
In the CLI provided by the UNIX or Unix-like OS, including macOS and different Linux distributions, the `ssh` command is usually installed from the beginning, which you can use to access another machine via SSH protocol. The way of using it is simple. You can open the CLI at your computer (e.g. the app Terminal in macOS), and then run the `ssh` command in the following manner:

```console
ssh <username>@<hostname>
```

Here you would need to make sure you have a username or account name available at the server, so as the name of the server (hostname). Also you need to make sure that your personally computer and the server is in the same network (e.g. to use the servers in ETHZ, you should make sure you are using the ETH network via WIFI or network cables in the campus, or you should use the ETH VPN), otherwise the computer won't be able to understand and connect to the right machine. Without other additional setup, the server will then respond by asking you for the password. If the server is working, your account is valid on the server, and the username and password you input match, you can likely see some welcoming message, and a new line with the prompt waiting for your command in the server.

For the students attending the Systems Genomics course, there is the student server for you to work on your drylab task. The server is called "bs-studentsvr04" as its hostname. The username to login the server is your ETH username, while the password being the one you use also for your email and most of the ETH services.

```console
ssh <username>@bs-studentsvr04
```

You would hopefully see the following, suggesting a successful login:
<p align="center">
<img src="img/ssh.png"/>
</p>

If you want to close the connection, type `exit`, or directly close the Terminal window.

#### If you are a Windows user
For the Win users, using `ssh` is unfortunately not as simple, as it is not natively provided in both the GUI and its CLI apps. Therefore, you need to firstly retrieve a SSH client app so that you can use it to connect to and access the server. There are quite some SSH client apps available. One of them, which is small but powerful and famous, is PuTTY. It is available here: https://www.putty.org/. From its download page, you can download the MSI ("Windows Installer") of the right system config (in most of the time should be 64-bit x86), and then execute the downloaded file when it is done. After finishing the installation, you should be able to find the PuTTY app from the Start menu.

The PuTTY app has a GUI as following:
<p align="center">
<img src="img/PuTTY.png"/>
</p>

The simplest and fastest way to start with is to put in the hostname of the server (bs-studentsvr04 for the student server for the Systems Genomics course), and then click on the "Open" button or type the Enter key directly, and it would open a terminal-like window where there is the text asking you to type in your username. Type in your username followed by Enter to input, it would then ask for your password. Enter your password, and if everything is alright, you will see things basically the same as the `ssh` example, suggesting a successful login.

>***NOTE***
>For some advanced Windows users, you may have heard of or been using the WSL (Windows Subsystem for Linux) system. Microsoft developed a compatibility layer for running Linux binary executables natively on a Windows system, and it is available in Windows 10 and Windows 11. From 2019, there is the MSL2 announced which introduced a real Linux kernel. It also supports different Linux distributions including Ubuntu, OpenSUSE and some others. If you have WSL in your Windows computer, you don't need to use PuTTY in principle, as you should directly have `ssh` installed as a part of your WSL. In that case, what you need to do is to simply open the WSL (e.g. search for Ubuntu in the start menu if that's the one you have), and then run `ssh` directly just as in macOS or Linux.

No matter in which way, once you login the server successfully, you are ready to play around with it a bit with what have been shown above, and then continue to the next section.

#### More information about the bs-studentsvr04 server


### 1-3. Install the required tools with the help from conda
Now you have access to the server, and hopefully also know something about how to use it via command line. Now it's time to set up all the tools needed for the following data preprocessing and analysis. Here I summarize some software which will be introduced and/or used later.

|Software|Link|Function|Compatible OS|
|--------|----|--------|-------------|
|SRA-Toolkit|https://github.com/ncbi/sra-tools/wiki|Retrieve data from SRA|UNIX/Unix-like, Win|
|SRA Run Selector|https://www.ncbi.nlm.nih.gov/Traces/study/|Interactive filter and selection of SRA entries to obtain their metadata and accessions|Online|
|FastQC|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|Quality control for the FASTQ files|UNIX/Unix-like, Win|
|Cutadapt|https://cutadapt.readthedocs.io/en/stable/index.html|Find and remove unwanted sequence from sequencing reads|UNIX/Unix-like, Win|
|STAR|https://github.com/alexdobin/STAR|RNA-seq read mapping|UNIX/Unix-like|
|kallisto|https://pachterlab.github.io/kallisto/|RNA-seq read pseudomapping|UNIX/Unix-like|
|Samtools|http://www.htslib.org/|View and manipulate SAM/BAM files|UNIX/Unix-like|
|RSEM|https://deweylab.github.io/RSEM/|Expression quantification|UNIX/Unix-like|
|R|https://www.r-project.org/|Commonly used programming language and analytical framwork for statistics|UNIX/Unix-like, Win|
|DESeq2|https://bioconductor.org/packages/release/bioc/html/DESeq2.html|Differential expression analysis|R package|

The SRA Run Selector is a webtool that you can simply access using your browser (e.g. Google Chrome), so no installation is needed and you can simply use your personal computer for it.

Some of those tools do not require too much effort to set up. They are either implemented with a cross-platform programming language (e.g. Java), or pre-compiled using the same or compatible system as the system used in the server. For those software, you can directly download it from the website to somewhere in the server (e.g. with `wget`), decompress it if needed (e.g. with `unzip` for .zip files, `gzip` for .gz files, `tar` for .tar, .tar.gz and .tar.bz files), and then there is the executable file available to run. This category includes FastQC (Java-based) and SRA-Toolkit (pre-compiled).

```console
cd [the students folder]
mkdir tools
cd tools
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastq_v0.11.9.zip
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64.tar.gz
tar xf sratoolkit.3.0.0-centos_linux64.tar.gz
wget http://opengene.org/fastp/fastp.0.23.1
```

For the others, the installation from scratch would not be the most pleasant work in the world. Unlike in Windows or macOS where one usually just need to run the installor of a software and follows the step-by-step guide to have it installed, a more common situation when using a Linux server as an ordinary user is that if you want to install a software, you need to download the source code, compile it, install the resulted executable command to somewhere that you have the access permission, and then tell the system to include that place into the searching path for command. Even when you don't run into any problem of dependency (e.g. one software may need another one which is not available in the system, so you would have to install the other one first), this is also a pain in the arse. If you are a system admin of the server (i.e. "sudoer", meaning the users who can use the command `sudo` to act as a root user that in principle can do anything on the system), you can probably install them quite easily by using the package manager repositories of the Linux system (if they are available). However, that would be unlike the case, and for sure won't happen if you use the student server for the Systems Genomics course. And even if you have the permission, it is usually the guideline to only install a software or package using the sudo permission when there is no alternative, so that the system can stay at the minimal scale for security and robustness.

Luckily, now we have conda.

[Conda](https://docs.conda.io/en/latest/) is an open source package management system and environment management system that runs on most of the OS including Windows, macOS and Linux. It has huge software repositories (called channels), and relying on them it quickly installs, runs and updates packages and their dependencies. When being asked to install a package or software, conda checks its/their availability in the channels, and if they are available, it retrieves the dependencies of the requested software to make sure of their availabilities in the environment; if any of them being missing, it is added to the list of things to install. The programs available in those channels are precompiled, so no compilation is needed to happen locally, which saves a lot of time and effort as well.

To use conda, we need to firstly install conda in the server. We can download miniconda, a free minimal installer for conda, to the server, and then run it. More information about miniconda can be found here: https://docs.conda.io/en/latest/miniconda.html.

```console
cd [the students folder]
cd tools
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh
bash Miniconda3-py37_4.12.0-Linux-x86_64.sh
```

And then simply follow the prompts on the screen. Most of the settings you can simply leave it default without any change. However, if you are using the bs-studentsvr04 for the Systems Genomics course, you should change the default installation location (your home folder) to your scratch folder. This also applies to many other computing servers, where limitations are set to how much data and/or file numbers you can store in the home folder.

Once the installation is done, you shall quit the login session to the server, and then log in it again. Afterwards, you can check whether the conda is successfully installed and set up by simply checking where your Python interpreter is.

```console
which python
```

As Python is included in miniconda, you should see the Python interpreter locates at your scratch folder where you install conda to. On the other hand, if what you see is the Python preinstalled in the Linux system, e.g.

```console
$ which python
/usr/bin/python
```

your conda is either not installed successfully, or not yet properly initialized. If you are sure that you see no error message during the installation process, you can try to initialize the conda set up again by `conda init`. Afterwards, quit your login session and start a new one, and then check again. If it still fails, you shall try to install conda again, or ask for help.

And once you have the conda installed and set up, it would be just one command to have the remaining tools installed.

```console
conda install -c bioconda -c conda-forge cutadapt star kallisto samtools rsem
```

Indeed, you can install FastQC, SRA-Toolkit and FASTX-Toolkit also with conda
```console
conda install -c bioconda -c conda-forge fastqc sra-tools fastp
```

It would ask you to confirm the installation of not only the four requested software but also all the dependencies. Once everything is finished, you can use the `which` command to make sure those tools are installed (e.g. `which STAR`).

>**NOTE**
>Conda has more functionalities that just installing software more easily. It is an environment manager, meaning that you can create, manage, delete and use different environment where different software are installed while making sure those different environments would not affect each other. This is especially useful when you want to use certain versions of software in some scenarios but not the others. Having too many software installed in the same environment also makes it more difficult and time-consumping to resolve the dependencies when installing new packages or upgrading the existing ones. Therefore, to have different environments set up for different purposes could be a very useful guideline in the future.

### 1-4. Get the public RNA-seq data from SRA
Now we have the computational environment ready for the preprocessing. We just need to data to start. As one can easily expect, there are two sources of data: the in-house data that are generated freshly by the lab for specific biological questions, and the public data which have been released and used for answer certain questions, but can be reanalyzed solely or together with other data for the same or related questions.

There are several huge repositories in the world for high-throughput sequencing data. The major players include [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) by NCBI in the US, [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) by EMBL-EBI in the UK, and [DDBJ Sequence Read Archive (DRA)](https://www.ddbj.nig.ac.jp/dra/index-e.html) by the DDBJ Center in Japan. These three repositories are also members of International Nucleotide Sequence Database Collaboration (INSDC), and cross-backup each other, meaning that data submitted to any one of the three databases are also accessible from the other two.

>**NOTE**
>While majority of the high-throughput sequencing data are archived in these three databases, there are also other emerging sequencing data repositories, though most of them are regional, mostly used by researchers in the host country. Examples include [Genome Sequencing Archive (GSA)](https://ngdc.cncb.ac.cn/gsa/) by NGDC, and [CNGA Sequence Archive (CNSA)](https://db.cngb.org/cnsa/) by CNGB, both located in China.

In this tutorial, we will retrieve the data we need from SRA. This is not only because the data we are going to use are archived at SRA, but also because of SRA-Toolkit which provides a simple way to download data given the accession numbers of the needed samples.

#### The example data set
The example data set used in this tutorial is from the paper [Comprehensive transcriptome analysis of neocortical layers in humans, chimpanzees and macaques](https://www.nature.com/articles/nn.4548) published in *Nature Neuroscience* in 2017. In the paper, the authors generated RNA-seq data representing the layer structure of prefrontal cortex in human, chimpanzee and rhesus macaque brains, aiming to identify human-specific transcriptome changes, which implies potential organizational rearrangement that specifically happened to human during evolution, and therefore contributed to the human specific features.

As this tutorial is not to reproduce the result presented in the paper, but to introduce the general procedure of preprocessing and analyzing RNA-seq data, we are not going to use the whole dataset, but only the subset of human samples representing the purely sampled layers.

To get the data, we firstly look at the paper. Nowadays most of the papers that involves large scale sequencing data also publish the accession number of the raw data. In most of the time, this information is included in sections like "Data availability", "Data and code availability" or "Accession codes". In this paper, there is indeed the "Accession codes" section:

```
Accessions
BioProject
PRJNA299472

Sequence Read Archive
SRP065273
```

#### Select the wanted samples via SRA Run Selector
From there we know that the data was deposited in SRA, with the accession number SRP065273. Of course, this includes all the data mentioned in the paper, and therefore not only the human samples, but also the other species. In addition, the paper presented two data sets, with one being the main one covering all the three species (DS1) and a second one with human and rhesus macaque only for verification (DS2), and both data sets are involved in the same accession. How shall we easily get the data subset that we want now, i.e. the DS1 human samples that represent purely sampled layers?

We can start with using [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/). It is an online tool by NCBI which can search for submitted SRA runs (usually means one round of sequencing of one sample) that are involved in the given accession. What's more it can do further interactive filtering based on certain metadata information of the samples.

First, we search for everything under the accession SRP065273.
<p align="center">
<img src="img/run_selector_1.png">
</p>

The SRA Run Selector outputs all the related sequencing runs, with the linked metadata also included in the table. It also has the metadata and summary table of the accession. If we look at the left hand side, there is the "Filter List" block, where we can choose what metadata information to do subsetting.

<p align="center">
<img src="img/run_selector_2.png">
</p>

Let's select "Data_Set" from the Filters List box, then select only "ds1" from the appeared box. Next we further select "Organism" from the Filters List box, the then "homo sapiens".

<p align="center">
<img src="img/run_selector_3.png">
</p>

Now the item (meaning sequencing run) number is reduced from 273 to 72.

However, because of the way of the sample collection, not all the samples purely represent one single layer of the cortex. Such information is unfortunately not included in the submitted metadata that we can directly see from the sample table. However, it is included in the paper, by integrating Figure 2 (for relationship between layers and the aligned sections) with the Supplementary Table 1 (for the aligned sections represented by each sample). Also since some layers are sampled more frequently than others due to the differences on their thickness,  Now we have a list of the 25 samples we will use.

|Individual|Sample|Layer|
|-----|-----|-----|
|DS1_H1|DS1_H1_01|L1|
|DS1_H1|DS1_H1_03|L2|
|DS1_H1|DS1_H1_06|L3|
|DS1_H1|DS1_H1_07|L4|
|DS1_H1|DS1_H1_10|L5|
|DS1_H1|DS1_H1_13|L6|
|DS1_H1|DS1_H1_18|WM|
|DS1_H2|DS1_H2_01|L1|
|DS1_H2|DS1_H2_02|L2|
|DS1_H2|DS1_H2_06|L3|
|DS1_H2|DS1_H2_08|L4|
|DS1_H2|DS1_H2_11|L5|
|DS1_H2|DS1_H2_14|L6|
|DS1_H2|DS1_H2_18|WM|
|DS1_H3|DS1_H3_02|L2|
|DS1_H3|DS1_H3_05|L3|
|DS1_H3|DS1_H3_06|L4|
|DS1_H3|DS1_H3_10|L5|
|DS1_H3|DS1_H3_13|L6|
|DS1_H3|DS1_H3_18|WM|
|DS1_H4|DS1_H4_03|L3|
|DS1_H4|DS1_H4_06|L4|
|DS1_H4|DS1_H4_09|L5|
|DS1_H4|DS1_H4_12|L6|
|DS1_H4|DS1_H4_18|WM|

>**NOTE**
>This is actually a bad example of metadata preparation when submitting data to repositories, although not the worst scenario when it is impossible to retrieve the critical information in any way (this is unfortunately quite common). As the data submitter of the paper, I truely appologize for this.

Now we can select those samples from the SRA Run Selector. After selecting only those samples, you can further switch on the "Selected" option at the "Select" block of the page to only show the selected samples in the item table and double check. Once confirmed, you can click on the "Metadata" and "Accession List" buttons at the "Download" column of the "Select" block, to download the complete metadata that you see in the item table (Metadata) and just the list of accessions of those samples (Accession List).

<p align="center">
<img src="img/run_selector_4.png">
</p>

#### Download the raw sequencing data in FASTQ format via SRA-Toolkit
Now it is time to download the data. As mentioned above, SRA-Toolkit provides the functionalities to download data from SRA, given the list of SRA sequencing run accessions which we just obtained. And of course, we would like to do the download directly on the machine that will be used for the following preprocessing (e.g. the bs-studentsvr04 server). Since you probably just downloaded the accession list to your personal computer, you need to upload it to the machine. There are different options for that.

You can just copy-paste the content of the accession list file. Open the accession list by any text editor on your computer, select all the context and then copy, and then login the server, go to your work folder (e.g. [the students folder]), and then create to paste all the content. For instance, you can use nano by `nano SRR_Acc_List.txt` to create the file "SRR_Acc_List.txt", then in the editor paste the text just copied with cmd+C for macOS or right-click at the PuTTY window if you are a Windows PuTTY user. After seeing all the accession numbers being pasted to nano, press ctrl+X (as indicated at the footnote menu, `^X Exit`. The character ^ in front of a Latin alphabet means pressing both the ctrl button and the other character button at the same time) for exit, and then press `Enter` to confirm the saving.

Alternatively, you can directly transfer the file to the server. This sounds overkilling in this case but it is good to know how to do it as you may need to copy something much bigger to the server. For most of the servers supporting SSH access, it also supports data transfer via SFTP (Secure FTP, or File Transfer Protocol), or the `scp` command mentioned above if your personal computer is on macOS or Linux. The way of using `scp` is very similar to `cp`. The only difference is that you need to include your username and the hostname of the remote machine in the file name so that the command knows that the file source or target is at a remote computer. For example,

```console
scp Downloads/SRR_Acc_List.txt hezhi@bs-studentsvr04:/mnt/users/hezhi
```
In this case, the command will then ask for the password for logging in the bs-studentsvr04 machine with the username hezhi (this is my username. Don't forget to change it to yours). Type in the matched password (the same password as when logging in via ssh), and then the transfer will start, and the new copy will be at the `/mnt/users/hezhi` directory in the bs-studentsvr04 server.

Once you get the accession list in the server, you can do the data download using SRA-Toolkit. The toolkit contains many different commands. Among them the most relevant ones include `prefetch`, `fastq-dump` and `fasterq-dump`. The `prefetch` commands can take a accession list file as the input and download the data of those accessions in the format of [SRA data format](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/) which the SRA repository uses to store sequencing data. It is however not the standard data format of sequencing data that any genomic data processing tool will use. One can then use the `fastq-dump` command to convert the SRA files to the standard data format FASTQ, given all the downloaded SRA files with the glob pattern. 

```console
cd [student folder]
mkdir rawdata
cd rawdata
prefetch --option-file ../SRR_Acc_List.txt
fastq-dump --gzip --split-3 SRR*/*.sra
```

>**NOTE**
>* The `prefetch` command saves the downloaded data of each sequencing run (i.e. every unique SRR accession) in a separated directory named by the SRR accession.
>* When using `fastq-dump`, it is important to add the `--split-3` option, so that when the sequencing data is paired-ended, the two mates are automatically split.
>* By default, the `fastq-dump` output the final FASTQ files in plain text which can be very large. Adding the `--gzip` option directly compresses the FASTQ files with gzip so that it becomes ~1/3 as big. It also means taking a bit longer time for the compression but to save the limited storage this is usually worthy.

The alternative way is to use `fastq-dump` or `fasterq-dump` to download and convert to FASTQ files directly without the SRA intermediates. These two commands can both do the same work but some differences also exist. For instance, `fasterq-dump` allows multiple downloading threads which can very likely speed up the download. `fasterq-dump` also set the `--split-3` option as default so that you don't need to worry about forgetting it. On the other hand, the `--gzip` option is not supported in `fasterq-dump`, so you would have to explicitly run `gzip` later on the downloaded FASTQ files. In general, it is recommended to use `fasterq-dump` than `fastq-dump` but it doesn't really matter too much. One issue for both the commands is that they don't support the accession list file as the input, but expect the explicit given accessions, meaning something like
```console
fasterq-dump --threads 5 --progress SRR2815952 SRR2815954
```

>**NOTE**
>In the `fasterq-dump` command, the `--threads` option specifies the number of downloading threads (by default 1, so no threading). The `--progress` option prints the downloading progress to the screen so that you know how much it has gone easily.

This is a bit annoying as we don't want to type in many SRR accession numbers one by one manually. The solution here is to use the piping feature mentioned above together with the `cat` command to print content of a given list, and the `xargs` command that convert standard input into arguments of another command
```console
cat ../SRR_Acc_List.txt | xargs fasterq-dump --threads 5 --progress
gzip *.fastq
```

Once the download is finished, you can list the files in your working directory and see whether you can all the files as expected. They should all be named as [SRR Accession].fastq.gz.

```console
$ ls -l *.fastq.gz
-rwxrwx--- 1 hezhi@d.ethz.ch bsse-treutlein@d.ethz.ch 102401267 Aug 29 15:17 SRR2815952.fastq.gz
-rwxrwx--- 1 hezhi@d.ethz.ch bsse-treutlein@d.ethz.ch 533150939 Aug 29 15:17 SRR2815954.fastq.gz

...[skip the other lines]...

$ ls -1 *.fastq.gz | wc -l
25
```

>**NOTE**
>This example data set is single-ended RNA-seq data, therefore each sequencing run has only one FASTQ file. For the paired-ended RNA-seq data, each sequencing run is split into two FASTQ files, each for one mate. The filenames would be {SRR Accession}_1.fastq.gz and {SRR Accession}_2.fastq.gz.

#### An introduction to the FASTQ and FASTA data formats
FASTQ is the standard data format used to store sequencing data. It was originally developed at the Wellcome Trust Sanger Institute to incorporate quality information into the FASTA format, which is the standard data format developed to represent nucleotide or amino acid sequences where every nucleotide or amino acid is represented by one single character (e.g. A/T/C/G for DNA sequences).

In a FASTA file, one or multiple sequences are represented, each in a unit of two components. The first component is strictly one line called "description line", which starts with the character ">". The description line gives a name and/or a unique identifier for the sequence, and potentially also additional information. The second component is the sequences in one or multiple lines, one character for one nucleotide or amino acid. In genomics, the most commonly usage of FASTA is to store the sequences of genome, annotated transcriptome and peptide products. The FASTA file looks like the following:
```
>sequence A
GGTAAGTCCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATAT
TCTGTTGCCAGAAAAAACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC
>sequence B
GGTAAGTGCTCTAGTACAAACACCCCCAATATTGTGATATAATTAAAATTATATTCATAT
TCTGTTGCCAGATTTTACACTTTTAGGCTATATTAGAGCCATCTTCTTTGAAGCGTTGTC
TATGCATCGATCGACGACTG
```

A FASTQ file looks very relevant but different from a FASTA file. It is used to describe one or multiple sequences, but in most of the time multiple (usually tremendous amount of) sequences. Every sequence is represented in a unit of four lines:
* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.

It looks like the following:
```
@SRR2815952.1 1 length=100
AGACGAGACCTACTGCATTGATAACGAAGCTCTCTACGACATTTGCTTCAGAACCCTAAAGCTGACCACGCCCACCTATGGTGACCTGAACCACCTGGTG
+SRR2815952.1 1 length=100
@?@B?@DDDDHDCGHGHIIEHGIIFEHI@?BFGABDFGHAGGIIIHIIIGGIIIBHIIFGIHACHHEEBEBCCDD@ACC:>>@CDDD?CCD(<<?A?A@C
@SRR2815952.2 2 length=100
TGGGGTTTCACCATGTTGGCCGGGCTGGCCTCGAACTCCTGACCTTGTGATGCACCCACCTCGGCCTCCCAAAGTGCTGGGATTTACAGGCGTAAGCCAC
+SRR2815952.2 2 length=100
CCCFFDFFHHHHHJJJJJJJJHGIJJJJJJJJJIJJIJIIGIJJJJJIIJJJJJJJHHHFFFFDDDDDDDDDBC@>@DDDDDDDDDEDCCBDDDDDDDDD
```

FASTQ is a great format to represent short sequences with quality values, and therefore high-throughput sequencing data where the quality values can represent the sequencing quality. In the Illumina sequencing, a Q score is used to represent the quality of each sequenced nucleotide. The Q score is defined as
$$Q = -10 \log_{10} e$$
Here, $e$ is the estimated probability of the base call being wrong. Therefore, a higher Q score means lower probability of sequencing error (or more precisely base calling error), and therefore higher quality. A Q score of 20 (Q20) means an error rate of 1 in 100, or accuracy of 99%. A Q score of 30 (Q30) means accuracy of 99.9%. At such a level, the base is considered to be accurate enough that all bases of a read is likely to be correctly called if all the bases in a read reach such a level of quality. Currently Q30 is considered a benchmark for quality in next-generation sequencing (NGS).

To represent the Q scores of all the based in a sequence in the FASTQ format, the numeric Q score of each base is encoded into a compact form based on the [ASCII codes](https://www.ascii-code.com/). Basically, the estimated Q score represented as an integer is represented by the character with the ASCII code equal to $Q+33$. It has the 33 component because Q is strictly non-negative with the minimum of 0 meaning 0% accuracy, and 33 is the smallest ASCII code that defines a one-character symbol (!). 

## Preprocessing of RNA-seq data
Now we have the tools ready, and the data ready. It is time to move on to the next step, to preprocess the RNA-seq data. In general it contains the following steps:
1. Quality control
2. Read mapping or pseudomapping
3. Gene expression quantification for samples
4. Generate the expression matrix with sample metadata table for the following analysis

### 2-1 Quality control of RNA-seq data
Before actually processing the data, it is important to make sure that the data is of high quality.

The quality of a RNA-seq data contains different aspects. The sequencing quality, represented by the quality score for each base in each read, is one of the very first thing one should consider. It varies from one base to another base, and from one read to another. We can therefore look at the distribution of the Q score per loci across all reads. Usually the sequencing reads have relatively low quality on the two sides (especially the end) and higher quality in the middle, and this is the reason why we look at the per loci quality distribution across reads. Another reason to look at that is because if it is indeed the case that the start and/or the end of reads systematically have too low quality, one can easily fix it by trimming the first and the last bases of all the reads. This will be mentioned in more details later.

Other quality metrics are more related to the sample and cDNA library quality. For instance, the adapter ligantion is a critical part during the cDNA library preparation and cases can happen that multiple adapters are ligated to the sequence and therefore become parts of the sequenced reads. This would introduce troubles later when we need to locate the transcript(s) that the read represents, as the extra adapter sequences are not a part of the transcripts and would introduce a large mismatch. Another example is the complexity of the data. Ribosomal RNA (rRNA) makes up about 80% of cellular RNA, while the rRNA genes only makes up 0.5% of the human genome, and their abundances are not very relevant to many biological processes to study. Therefore, there is usually a step of mRNA enrichment (by Oligo-T sequences) or rRNA depletion to effectively capture the more informative non-rRNA transcript fractions. However, this is not always working efficiently and the cDNA library would present a low complexity if the rRNA portion is not effectively reduced. Also, there could be reads representing pure ligation products of multiple sequencing adapter, which also dramatically reduce the library complexity.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a tool providing a simple way to do some quality control checks on the sequencing data. It checks different aspect of data quality and provides a graphical report so that one can intuitively get the idea about the data quality. The tool contains both GUI and CLI, therefore one can choose to either run it by clicking buttons with mouse or with commands.

To run it in the command line, it is very simple. This is an example:
```console
cd [student folder]
cd rawdata
mkdir fastqc
fastqc -o fastqc *.fastq.gz
```
>**NOTE**
>The `-o` option in the `fastqc` command specifies the output directory to put the output files. By default it stores the output files to the same folder as the FASTQ files. In the example script, a new folder is created to store just the FastQC reports. This is not a must, but will likely help to make it more tidy.

If you prefer GUI, you can also directly run `fastqc` to open the GUI in the server. If everything goes well you would be able to see the FastQC window poping up. Then you can open one or more FASTQ files by selecting from the menu `File` -> `Open...`. After getting the report that's shown in the window, you can save it, also two files just as the command line version, by choosing from the menu `File` -> `Save report...`.

>**NOTE**
>Please note that to see the GUI at your personal computer, you need to make the X11 forwarding possible. [X11](https://en.wikipedia.org/wiki/X_Window_System) is a windowing system for bitmap displays, proving the basic framework for a GUI environment and is commonly used in Linux. To allow the X11 to work via SSH-based remote access, there are several things one would need to do:
>* For macOS users, make sure to add the `-Y` option when running the `ssh` command. It means to use `ssh -Y <username>@bs-studentsvr04` to login the bs-studentsvr04 server, for instance. Also, you need to make sure that XQuartz, the X11 system that runs on macOS, is installed.
>* For Windows users using PuTTY, it is more complicated. First of all, you need to make sure that [Xming](http://www.straightrunning.com/XmingNotes/) or any other X11 server for Windows is installed, and being opened. And then in PuTTY, after filling in the server information at the starting page, from the option menu on the left hand side, choose "Connection" -> "SSH" -> "X11". There you shall click the "Enable X11 forwarding". Afterwards click the "Open" button to start the connection.
>* The newest Xming is not free, but needs £10 "donation" to get the download password. Its older version, however, is available at [sourceforge.net](https://sourceforge.net/projects/xming/files/Xming/6.9.0.31/Xming-6-9-0-31-setup.exe/download) for free.

For each FASTQ file there are two output files generated by default. Assume the FASTQ file is called `[filename].fastq` or `[filename].fastq.gz`, then one output file is called `[filename]_fastqc.html` and the other one called `[filename]_fastqc.zip`. The HTML file is a report which can be opened with your browser, and it contains the summary of all the quality metrics, as well as a grade by the software on each perspective whether it is passed, warning, or failed. The ZIP file, once decompressed, contains also the HTML report, as well as plain text files with the same information.

This is an example screenshot of the HTML report when being opened in the browser:

<p align="center"><img src="img/fastqc.png" /></p>

It is important to mention that some of the grades are made assuming the data to be whole-genome DNA sequencing. For instance, the "Per sequence GC content" compares the distribution of G/C bases proportion per read to a theoretical distribution derived from the whole genome, which is expected to be different from the GC content of transcriptome. From my personal experience, the "Per base sequence content" and "Per sequence GC content" are the two sections that easily get the warning for failed grade for RNA-seq data, but can be ignored if other sections are fine. In addition, the "Sequence Duplication Levels" is another section that could give out warning of RNA-seq data, while it may or may not be a problem that needs to be solved later.

Meanwhile, the sections that I would suggest to pay attention to for RNA-seq data include "Per base sequence quality", "Sequence Duplication Levels", "Overrepresented sequences" and "Adapter Content". They represent potential We will need to try to fix the problem if they get a failed grade:
* If any read locus shows low quality (e.g. median <20 or even <10) in the "Per base sequence quality" section, especially at the two ends, we should try to trim them if the low-quality part is large (>10 bases), either by all reads removing the same number of bases or different number per read based on the quality scores.
* Since different transcripts have very different abundance, to make sure that very lowly expressed transcripts are also detected, it is possible that the highly expressed transcripts are over-amplified and/or over-sequenced, resulting in warning of "Sequence Duplication Levels". In this case, a de-duplication step may be wanted to collapse the identical reads into one.
* For standard mRNA-seq data with oligoT enrichment, problems of "Overrepresented sequences" and "Adapter Content" often come together and represent the adapter ligation issue mentioned above. We can try to cut the adapter sequences from reads later.

For the example shown in the screenshot above, we don't need to do anything as it looks all good.

**IMPORTANT NOTE**: It is not always necessary to do anything here even if problems were found, especially those related to base quality. For instance, many up-to-date software being used later for read mapping (e.g. STAR) has implemented a soft trimming mechanism to deal with low-quality bases at the end of a read.

#### (Optional) Trimming and deduplication
If it is really needed, fixing the first and the third problems requires cutting off parts of some/all reads in the data. There are quite some tools providing the functionality to do that, including [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) and [fastp](https://github.com/OpenGene/fastp). The fastp tool, although less famous and commonly used than Cutadapt, provides the de-duplication function altogether so that one can do trimming+deduplication at the same time.

Let's firstly look at Cutadapt. It has a lot of functionalities related to do all kinds of trimming in pretty complicated manners. You can get all the details from its [online user guide](https://cutadapt.readthedocs.io/en/stable/guide.html). As an easy example, assuming the data being single-ended (SE), we can trim the adapter sequences in the data by using the following command:
```console
cutadapt --adapter=AGATCGGAAGAG --minimum-length=25 -o SRR2815952_trimmed.fastq.gz SRR2815952.fastq.gz 
```

For paired-ended (PE) reads, it is very similar 
```console
cutadapt -a=AGATCGGAAGAG -A=AGATCGGAAGAG --minimum-length=25 -o read_1_trimmed.fastq.gz -p read_2_trimmed.fastq.gz read_1.fastq.gz read_2.fastq.gz
```

>**NOTE**
>* In both examples, the sequence AGATCGGAAGAG is the first 12 bases of the Illumina Universal Adapter. It is exactly one of the sequences that FastQC uses to evaluate adapter content. If there is a different adapter sequence being used, the related options should be changed accordingly.
>* If the data is PE, make sure to provide both reads to trigger the PE trimming mode, so that the two paired reads are kept or removed together in a pair. This is important. For PE data the two reads are stored in two separated FASTQ files, but it is strictly required that reads in the two files have the same number of reads which are paired in the same order.

As you can see above, the `cutadapt` command expects only one sample per time. When you have many samples, you need to either run it one-by-one manually, or you can rely on for-loop in the Bash scripting. For example, this is the script to apply the same trimming to all the fastq.gz files in the current folder, and store the trimmed FASTQ files in the new subdirectory called `trimmed`:
```console
mkdir trimmed
for file in *.fastq.gz; do
  cutadapt --adapter=AGATCGGAAGAG --minimum-length=25 -o trimmed/$file $file
done
```

Similar process can be done by using fastp, with the deduplication operation also applied.
```console
fastp --dedup --adapter_sequence AGATCGGAAGAG -l 25 -i SRR2815952.fastq.gz -o SRR2815952_trimmed_dedup.fastq.gz
```

And of course we can wrap it up and apply to all the .fastq.gz files in the folder using the for-loop in Bash
```console
mkdir trimmed_dedup
for file in *.fastq.gz; do
  fastp --dedup --adapter_sequence AGATCGGAAGAG -l 25 -i $file -o trimmed_dedup/$file
done
```

Also to keep in mind that fastp is able to do more complicated manipulations and examples are shown in its [github page](https://github.com/OpenGene/fastp). It also provides a QC summary, not as comprehensive as FastQC does but still reasonable. However, we won't go into those details in this tutorial.

## 2-2 Read mapping/pseudomapping and quantification
### 2-2-1 Read mapping via STAR and data quantification
Once the quality of the data is confirmed, we need to convert those millions of reads per sample into the gene- or transcript-level quantification. This would need the assignment of reads to genes or transcripts. To do this, the mostly common first step is for each read, to look for the genomic region that match with the read, given the complete genomic sequences. The identified region is then most likely the region being transcribed and generate the sequenced read in the end. This step of looking for the matched genomic regions for reads is called read genome mapping or alignment.

There are different tools, or aligners, that have been developed for this purpose. The most famous examples include [Tophat/Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)/[HISAT2](https://daehwankimlab.github.io/hisat2/) and [STAR](https://github.com/alexdobin/STAR). As the commonly used modern aligners, HISAT2 and STAR shares quite some features, such as their high-efficiency, and their support of soft-trimming for low-quality bases at the ends of reads. They also have their own adventage and disadventage. HISAT2 uses fewer computational resource than STAR (particularly memory) and has better support for SNPs (single-nucleotide polymorphism) that in the same locus on the genome different individuals can have different nucleotides. On the other hand, STAR is suggested to provide more accurate alignment results. It also supports varied ways for the next step to quantify transcript abundance. In this tutorial, we will use STAR to map the FASTQ files we retreived from SRA to the human genome.

#### Brief introduction to the STAR aligner
Before STAR was developed and got widely acknowledged, there had been other aligners being developed, with the most commonly used example being Tophat by Cole Trapnell when he was in his PhD in University of Maryland (he is now an Associate Professor in University of Washington). Those tools were great and used by many studies using RNA-seq whichhad shown its great potential but was not yet fully mature as it is today. The main problem of those tools before STAR was their speed. They might be good enough when there were several or dozens of samples to process, but not for the research project with huge consotia effort such as ENCODE (Encyclopedia of DNA Elements), which generated RNA-seq data for hundreds or even thousands of samples.

To solve the speed issue was one of the major motivations that STAR (Spliced Transcripts Alignment to Reference) was developed in 2009 by Alexander Dobin in Cold Spring Harbor Laboratory (he is now an Assistant Professor in CSHL). The [STAR paper](https://academic.oup.com/bioinformatics/article/29/1/15/272537) was published in 2013 in Bioinformatics. Until now, the tool is still under active improvement and maintenance.

In brief, STAR uses a two-step procedure to achieve the high-speed alignment. The first step is seed search. For each read, STAR firstly searches for its longest sub-sequence that perfectly matches at least one locus on the reference genome. This fragment is called maximal mappable prefixes (MMP). After getting the MMP for the read, STAR searches MMP again but only for the unmapped portion (i.e. the parts outside of the MMP) of the read. These two MMPs obtained by the sequential search are also altogether called *seeds*, and that's why this step is named "seed search". This sequential search not only make it straightforward to deal with splice junctions where a read contains sequences from two separated exons in the genome, but also greatly speed up the alignment as searching for the entire read sequence is no longer needed.

<p align="center">
<img src="img/star_seed_search.jpeg" /><br/>
<sub><i>Figure 1 in the STAR paper</i></sub>
</p>

To further deal with the possible mismatches (due to sequencing errors on reads, SNPs, point mutations, errors in the reference genome, etc.), when the MMP search doesn't reach the end of the read, the MMPs will serve as anchors in the genome that can be extended to allow for alignments with mismatches. If the extension procedure does not yield a good genomic alignment (due to poor quality at the ends of reads, poly-A tails, adapter sequence ligation, etc.), a soft trimming is applied (the remaining read sequence is ignored although not physically cut off).

After the seed search is done, STAR applies the second step, which is clustering, stitching and scoring. Seeds are firstly clustered based on proximity to a set of ‘anchor’ seeds. Then, seeds that map close enough around the anchors (so that it can be still considered to be an intron) are stitched together. In this way, different seeds of a read which are from different exons can be stitched. The stitching is guided by a local alignment scoring scheme to penaltze mismatches, insertions, deletions and splice junction gaps. The stitched combination with the highest score is chosen as the best alignment of a read.

More details information are available in the STAR paper (technical details in its Supplementary Materials).

#### Download genome sequences and create indexed genome for STAR
The seed search in STAR is implemented through uncompressed suffix arrays (SAs). We won't go into the details how SAs works, but STAR needs the reference genome to be represented in the form of SAs so that it can apply its seed search for the reads. Therefore, the first step of using STAR is actually to have your reference genome ready, and then use the "genomeGenerate" function in STAR to prepare SAs for it.

So we need to firstly retrieve the reference genome. Specifically for the example data set, we need a human reference genome.

There are different places where we can obtain the genome sequences of human and other species. The [UCSC Genome Browser](https://genome.ucsc.edu/) provides probably the most and up-to-date resources of genome sequences with selected gene annotations, which can be both browsed via the online browser and downloaded. It of course includes the human genome, which can be downloaded via "Downloads" -> "Genome Data" at the header menu bar.

<p><img src="img/ucsc_genome_browser.png" /></p>

In the download page, data are grouped by species in mostly alphabetical order, except human as the first one. The most popular species also have links at the very top so that one doesn't always need to scroll down that much. At the [human section](https://hgdownload.soe.ucsc.edu/downloads.html#human), you can see the links of different human genome data, which are further grouped by different human reference genome versions. Here we want to the newest human reference genome (GRCh38/hg38), and from that section we can choose "Genome sequence files and select annotations" and then "Standard genome sequence files and select annotations". That leads to a page with introduction and file description on top and links to different files for download if you scroll down to the bottom. Among all the provided files, what we need right now is [hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), which is the gzip-compressed FASTA file containing sequences of all the human chromosomes (as well as the scaffolds and contigs that haven't yet been integrated into the chromosomes).

So now you can choose to download the file directly to your computer, and then transfer it to the server using the methods mentioned above (`scp` or SFTP), or directly download the data to the server. For the latter, you can use the `wget` command:
```console
cd [student folder]
mkdir genome
cd genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

STAR expects decompressed FASTA file(s) for the reference genome, so you shall decompress it before moving on:
```console
gzip -d hg38.fa.gz
```

Now we can build the genome index of the human reference genome for STAR
```console
mkdir star-index
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles hg38.fa
```

This command tells STAR to use 10 cores in the server to build the genome index for the genome in the FASTA file hg38.fa, and then store the resulted indexed genome in the newly created star-index directory. With this setting, it takes about one hour for the human genome to finish the indexing, with about 30GB RAM needed.
```
$ STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles hg38.fa
        STAR --runThreadN 10 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles hg38.fa
        STAR version: 2.7.10a   compiled: 2022-01-14T18:50:00-05:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
Sep 01 12:16:25 ..... started STAR run
Sep 01 12:16:25 ... starting to generate Genome files
Sep 01 12:17:25 ... starting to sort Suffix Array. This may take a long time...
Sep 01 12:17:41 ... sorting Suffix Array chunks and saving them to disk...
Sep 01 12:57:58 ... loading chunks from disk, packing SA...
Sep 01 12:59:26 ... finished generating suffix array
Sep 01 12:59:26 ... generating Suffix Array index
Sep 01 13:07:15 ... completed Suffix Array index
Sep 01 13:07:16 ... writing Genome to disk ...
Sep 01 13:07:21 ... writing Suffix Array to disk ...
Sep 01 13:08:08 ... writing SAindex to disk
Sep 01 13:08:14 ..... finished successfully
```

#### Mapping with STAR
Once the genome indexing is done, you are ready to map the reads to the reference genome with STAR. STAR has pretty good default mapping-related parameter settings, therefore, they can be kept as default for most ordinary RNA-seq data sets. Parameters that always need to be set properly are those related to the input and output files. This is the basic example script to run STAR on one of the sample in the example data set (no hurry to run, there is more to come :wink:):
```console
cd [student folder]
mkdir mapping
mkdir mapping/SRR2815952
STAR --genomeDir genome/star-index \
     --runThreadN 10 \
     --readFilesIn rawdata/SRR2815952.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix mapping/SRR2815952/
```

>**NOTE**
>The `\` character at the end of a line tells Bash that the next line is also a part of the same command. Without `\`, the end of a line automatically means the end of a command. This is a good way to do if you want to keep the script pretty when there are many options and/or arguments in a command, but it is not required. Essentially it has no difference at all to the version with everything in one line

Here are some explanations to the parameters:
* `--genomeDir genome/star-index`: specify the location of the indexed genome directory, i.e. the output of the previous step
* `--runThreadN 10`: use 10 cores for the alignment
* `--readFilesIn rawdata/SRR2815952.fastq.gz`: specify the input FASTQ file(s). If the data is PE with the two mates in data_1.fastq.gz and data_2.fastq.gz, put them both here separated by a space (`--readFilesIn data_1.fastq.gz data_2.fastq.gz`)
* `--readFilesCommand zcat`: specify the command used to display content of the input file (the `--readFilesIn` option). `zcat` is the command to decompress a .gz file with `gzip` and then print to the screen (stdout). `zcat file.gz` is equivalent to `gzip -cd file.gz`. Since the input file is gzipped, `zcat` needs to be specified as the way to print it. Otherwise, the FASTQ file would need to be decompressed beforehand
* `--outSAMtype BAM SortedByCoordinate`: specify the type of output file. The standard output is in SAM format, and BAM is the binary version of SAM format which greatly reduces the storage usage, while more time is needed for the format conversion. A BAM file can be unsorted (`Unsorted`), or sorted by coordinate (`SortedByCoordinate`).
* `--outFileNamePrefix mapping/SRR2815952/`: specify how the output files are named. What is specified here is the prefix of the output file. For instance, the default BAM file name is `Aligned.sortedByCoord.out.bam`. In the example, the output BAM file will be called `mapping/SRR2815952/Aligned.sortedByCoord.out.bam`. Of course, in this case it becomes necessary that the directory "mapping" exists in the current folder, and there is also the directory "SRR2815952" in the "mapping" folder. This is why the two `mkdir` commands are used before the `STAR` command

#### Quantification of gene expression
For RNA-seq data, read mapping is not what we ultimately need. In most of the time, what we actually need is the assessment of expression levels of different genes, and mapping is just the intermediate step. So how shall we then convert the read mapping results to gene expression values?

The most straightforward way, is to count the number of reads that are aligned to the exonic region of each gene. The more reads are aligned to the gene, the higher expression the gene has. Obviously, such raw read count values have a critical problem, that different RNA-seq data can have huge difference on sequencing depths. For instance, we have one sample with 100 reads aligned to a gene, while the same gene got 200 reads in another sample. Does it mean the gene has higher expression in the second sample? Not necessarily, as we might have 1 million reads in total for the first sample, while 10 million reads for the second sample. With 10-fold higher coverage for the second sample, we would also expect to see 10-fold as many reads mapped to a gene with the same expression level in the two samples; while in this case, the gene only has twice as many reads aligned to the gene. This suggests that this gene probably has much lower expression in the second sample instead.

What just mentioned is one issue when comparing the same gene across different samples. There is also an issue when comparing different genes even in the same sample. The standard bulk RNA-seq experiments usually include a random fragmentation step so that different parts of one long transcript can all the sequenced. On the other hand, it means that longer transcripts are more likely to generate more fragments, and therefore more reads. Therefore, two genes with the same number of reads aligning to don't necessarily mean they have similar expression levels, if their transcript lengths vary a lot.

To take into account those biases, people introduce different so-called normalization methods, to try to derive some metrics of expression levels which are more comparable among genes and samples. The simple but commonly used options include RPKM (Reads Per Kilobase Million reads) or FPKM (Fragments Per Kilobase Million reads) which are very similar with each other, and TPM (Transcript Per Million reads). The two methods, RPKM/FPKM and TPM are very similar, which is to calculate a scaling factor per sample per gene to correct the biases due to differences at sequencing coverage and gene lengths. The difference between them is the RPKM/FPKM considers the total number of reads in a sample as the proxy of sequencing coverage. Essentially it assumes different samples contain the same amount of nucleotides. On the other hand, TPM considers the total number of detected transcripts as the proxy of sequencing coverage, which means it assumes different samples share the same number of transcripts. Practically speaking, both methods firstly do a gene-length correction by dividing the read number at a gene by the gene length (or times a fixed scaling factor, e.g. $10^6$ afterwards). After that, RPKM/FPKM divides the resulted values by the total number of reads/fragments in the sample, while TPM divides the results by the sum of the scaled values across all genes.

For both RPKM/FPKM and TPM, there is still one critical issue: the gene length.


<style scoped> table { font-size: 0.8em; } </style>

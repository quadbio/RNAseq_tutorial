@settings{
table {
  font-size: 0.8em;
}
}

# Tutorial for bulk RNA-seq data preprocessing and analysis
#### Compiled by Zhisong He
#### Updated on 24 Aug 2022
### Table of Content
  * [Introduction](#introduction)
  * [Preparation](#preparation)
    * [1-1. Linux and Bash, your important buddies for RNA-seq data analysis](#1-1-linux-and-bash-your-important-buddies-for-rna-seq-data-analysis)
    * [1-2. Access the computing server](#1-2-access-the-computing-server)
    * [1-3. Install the required tools with the help from conda](#1-3-install-the-required-tools-with-the-help-from-conda)
    * [1-4. Get the public RNA-seq data from SRA](#1-4-get-the-public-rna-seq-data-from-sra)
  * Preprocessing of RNA-seq data
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

<span style="font-size:0.8em">*P.S. Unix is a family of multitasking, multiuser computer operating systems that derive from the original AT&T Unix, and are characterized by a modular design (Unix philosophy). Different components work together, with the kernel being the center to provide the most basic management and support. Unix is not free, but it inspired the development of many free Unix and Unix-like OS, and importantly, the GNU (meaning GNU's Not Unix) project and the Linux kernel which later on further derived into different Linux distributions (like different branches or versions of Linux). Among them, Red Hat Enterprise Linux (RHEL), Fedora, CentOS, SUSE Linux Enterprise, openSUSE, Debian, Ubuntu, Linux Mint and Arch Linux are among the most famous ones.*</span>

<span style="font-size:0.8em">*P.S.2. Actually the Apple macOS (since version 10.5 with 10.7 Lion as the exception) is a UNIX 03-compliant OS certified by The Open Group, so although not a Unix-like (meaning look like Unix but is not been certified, e.g. Linux) but a UNIX. This is also the reason why macOS in prefered in relative to Windows if you have to use your desktop/laptop for the task.*</span>

From our side as the computer end users, we don't need to care too much about how different OS work at the low level. However, we need to interact with the OS, for instance, to ask it to open a software or whatever, and that varies a lot from one OS to another. Windows and macOS look very different and we know that. When switching from one to the other for the first time, most of us would need quite a long time to get used to everything. However, that usually wouldn't be too difficult as both OS have pretty and straightforward graphical design for you to interact with the system. This is called a graphic user interface (GUI). Linux also have GUI, and different Linux distributions have different ones. However, what is different from Windows and macOS is that the GUI is not essential for Linux. For many computers running on Linux, especially those for high-performance computing, the GUI components are not even installed. Instead, people using those computers rely on CLI, which is short for command-line interface.

<p align="center">
<img src="img/GUI_vs_CLI.jpg"/>
<br/><span style="font-size:0.8em"><i>Image adapted from https://www.onetekno.my.id/2021/12/perbedaan-gui-dan-cli.html.</i></span>
</p>

CLI has a much steeper learning curve than GUI, and that's exactly the reason why GUI was developed in response to the complain to CLI. So why are we still using CLI heavily?That was also a question I had at the beginning, but after using CLI for a while and finally getting into it, I realized at least several reasons.

1. <ins>**The fancy graphic view comes with cost**</ins>. Computing and displaying all the visual changes on the screen needs resources (the computing unit CPU and/or GPU, memory, and also the storage). When doing high-performance computing which requires a lot of resource, however, we would for sure want to maximize the resource we can use for the real computation we want to do rather than just to display a mouse cursor moving around. Don't underestimate how much this would have needed, especially keep in mind that many computing servers and HPCs are not just been used by just one person at a time, and one person may have many things running at the same time. With a GUI for every user at least would use tremendous amount of computer resource.
2. While the GUI can only do what the developers implemented explicitly, <ins>**the CLI is easy to program to work out something more complicated or repetitive**</ins>. For instance, renaming 1000 files with a similar manner (e.g. to rename everything from A\*\*\*\* to B\*\*\*\* while keeping the \*\*\*\* part the same) can be quite easily done in CLI using the basic rename command together with basic programming skills in several lines of commands, but you would probably need to find a tool specifically implemented for batch renaming to do that in GUI. One can also combine different steps which use different tools together quite easily in CLI, much easier than using GUI. Such kinds of operations are very common when dealing with large data (e.g. preprocessing a RNA-seq data set with 10 samples using the same pipeline).
3. <ins>**Developing GUI for software needs a lot of effort**</ins>. Designing and implementing a nice GUI for a program can take as much time as, if not more than, implementing the functional part of the program. For tools which are mostly used by CLI-familiar users, it doesn't make sense from the developer's perspective to implement the GUI. If all the tools running on the OS only use CLI, the GUI of the OS is not really that useful anymore.

Learning and getting used to CLI would be pretty tough at the beginning, but it is far from doing something impossible, and you would probably like it once you get familiar with the commands and learn a bit of the simple programming skills. It takes time and effort of course, but stay patient and believe in yourself that you can do it! 

Before going to the HOW-TO part, do keep in mind that there are different CLIs for different OS (yes, Windows and macOS both have CLI as well). And there are different CLIs implemented for even the same OS. For Windows, there are the Command Prompt which emulates many of the command line abilities available in MS-DOS (Microsoft Disk Operating System, the old OS with only CLI by Microsoft), and the new PowerShell which provides extended capacities of the Command Prompt. For UNIX and Unix-like OS including Linux and macOS, there are different types of Shell including C shell (csh), Bourne Again shell (Bash), and Z shell (zsh), an extension of Bash. In most of the time, Bash is the default CLI for a Linux regardless different distributions, while zsh is currently the default shell for macOS.

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
<span style="font-size:0.8em">*P.S. There are two types of output to the screen, one is called "stdout" (standard output) while the other one called "stderr" (standard error). There is also the "stdin" (standard input) which is usually the input from the keyboard. Formally they are called standard I/O streams, I/O means input/output, and streams here represent the flow of information. For the two standard output streams, stdout is usually the real output while stderr is usually for verbose, warning, or error message. The simple `>` only save stdout to the file. If you want to save the stderr, `2>` instead of `>` should be used, where "2" represents stderr. And by the way, as you may have guessed, "1" represent stdout, so you can also use `1>` which is actually the same as using `>` directly.*</span>

Besides, there is another important feature of Bash, the pipes, indicating by `|`. It is a bit similar to redirection, but instead of saving the output to a file, it directly use the output of one command (which would be printed to the screen if you just run that command directly) as the input of the next command. This could be extremely useful to combine multiple commands for some complicated operations without the need to generate any intermediates.

```console
<command1> | <command2>
```

For example, the `ls` command has the option `-1` to print one file per line. Meanwhile, the `wc` command has the option `-l` to only output the number of lines in the given file or the piping input. We can therefore combine them two using the pipes

```console
ls -1 | wc -l
```

In this combination, the output of `ls -1` becomes the input of `wc -l`, so the number of lines in the `ls -1` is printed. As the line number of `ls -1` is the same as the number of files in the current working directory, the final output actually tells you how many files there are in the current folder.

<span style="font-size:0.8em">*P.S. One can do multiple piping to build a pipeline. Of course, that would require that every command being used in the pipeline supports the use of stdin as the input and can output their results to stdout so that it can be piped into the next command.*</span>


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

No matter in which way, once you login the server successfully, you are ready to play around with it a bit with what have been shown above, and then continue to the next section.

#### More information about the bs-studentsvr04 server


### 1-3. Install the required tools with the help from conda
Now you have access to the server, and hopefully also know something about how to use it via command line. Now it's time to set up all the tools needed for the following data preprocessing and analysis. Here I summarize some software which will be introduced and/or used later.

|Software|Link|Function|Compatible OS|
|--------|----|--------|-------------|
|SRA-Toolkit|https://github.com/ncbi/sra-tools/wiki|Retrieve data from SRA|UNIX/Unix-like, Win|
|SRA Run Selector|https://www.ncbi.nlm.nih.gov/Traces/study/|Interactive filter and selection of SRA entries to obtain their metadata and accessions|Online|
|FastQC|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|Quality control for the FASTQ files|UNIX/Unix-like, Win|
|FASTX-Toolkit|http://hannonlab.cshl.edu/fastx\_toolkit/|FASTA/FASTQ file manipulation|UNIX/Unix-like|
|STAR|https://github.com/alexdobin/STAR|RNA-seq read mapping|UNIX/Unix-like|
|kallisto|https://pachterlab.github.io/kallisto/|RNA-seq read pseudomapping|UNIX/Unix-like|
|Samtools|http://www.htslib.org/|View and manipulate SAM/BAM files|UNIX/Unix-like|
|RSEM|https://deweylab.github.io/RSEM/|Expression quantification|UNIX/Unix-like|
|R|https://www.r-project.org/|Commonly used programming language and analytical framwork for statistics|UNIX/Unix-like, Win|
|DESeq2|https://bioconductor.org/packages/release/bioc/html/DESeq2.html|Differential expression analysis|R package|

The SRA Run Selector is a webtool that you can simply access using your browser (e.g. Google Chrome), so no installation is needed and you can simply use your personal computer for it.

Some of those tools do not require too much effort to set up. They are either implemented with a cross-platform programming language (e.g. Java), or pre-compiled using the same or compatible system as the system used in the server. For those software, you can directly download it from the website to somewhere in the server (e.g. with `wget`), decompress it if needed (e.g. with `unzip` for .zip files, `gzip` for .gz files, `tar` for .tar, .tar.gz and .tar.bz files), and then there is the executable file available to run. This category includes FastQC (Java-based), SRA-Toolkit (pre-compiled) and FASTX-Toolkit (pre-compiled).

```console
cd [the students folder]
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastq_v0.11.9.zip
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64.tar.gz
tar xf sratoolkit.3.0.0-centos_linux64.tar.gz
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
```

For the others, the installation from scratch would not be the most pleasant work in the world. Unlike in Windows or macOS where one usually just need to run the installor of a software and follows the step-by-step guide to have it installed, a more common situation when using a Linux server as an ordinary user is that if you want to install a software, you need to download the source code, compile it, install the resulted executable command to somewhere that you have the access permission, and then tell the system to include that place into the searching path for command. Even when you don't run into any problem of dependency (e.g. one software may need another one which is not available in the system, so you would have to install the other one first), this is also a pain in the arse. If you are a system admin of the server (i.e. "sudoer", meaning the users who can use the command `sudo` to act as a root user that in principle can do anything on the system), you can probably install them quite easily by using the package manager repositories of the Linux system (if they are available). However, that would be unlike the case, and for sure won't happen if you use the student server for the Systems Genomics course. And even if you have the permission, it is usually the guideline to only install a software or package using the sudo permission when there is no alternative, so that the system can stay at the minimal scale for security and robustness.

Luckily, now we have conda.

[Conda](https://docs.conda.io/en/latest/) is an open source package management system and environment management system that runs on most of the OS including Windows, macOS and Linux. It has huge software repositories (called channels), and relying on them it quickly installs, runs and updates packages and their dependencies. When being asked to install a package or software, conda checks its/their availability in the channels, and if they are available, it retrieves the dependencies of the requested software to make sure of their availabilities in the environment; if any of them being missing, it is added to the list of things to install. The programs available in those channels are precompiled, so no compilation is needed to happen locally, which saves a lot of time and effort as well.

To use conda, we need to firstly install conda in the server. We can download miniconda, a free minimal installer for conda, to the server, and then run it. More information about miniconda can be found here: https://docs.conda.io/en/latest/miniconda.html.

```console
cd [the students folder]
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
conda install -c bioconda -c conda-forge star kallisto samtools rsem
```

Indeed, you can install FastQC, SRA-Toolkit and FASTX-Toolkit also with conda
```console
conda install -c bioconda -c conda-forge fastqc sra-tools fastx_toolkit
```

It would ask you to confirm the installation of not only the four requested software but also all the dependencies. Once everything is finished, you can use the `which` command to make sure those tools are installed (e.g. `which STAR`).

### 1-4. Get the public RNA-seq data from SRA


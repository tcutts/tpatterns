TPatterns should compile and install cleanly on any UNIX variant.  I have
tested it on:

* IRIX 6.4
* Linux
* Solaris 2.5
* MacOS Sonoma

# INSTRUCTIONS

1.  Install MPI and PCRE on your computer

1.  Type:
```
make
make install
```
3.  Create a suitable FASTLIBS file, and set the FASTLIBS environment variable to point to it.  A brief description of the format can be found in the manpage.
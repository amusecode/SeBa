**Installing & Compiling SeBa**

*Installing:*
- Download SeBa from github https://github.com/amusecode/SeBa

The only prerequisite software for SeBa is a C++ compiler.

*Compiling:*
To compile the code, several steps are needed.
Also after changing the code, I recomment repeating these steps.

On the commandline
Go to the SeBa directory, then do:
- make clean
- make
- cd dstar
- make

This will create an executable named SeBa in the subdirectory dstar
You can rename this executable to indicate the current version of SeBa.
You can move this executable to any directory.
When you run SeBa, an output file is created called ‘SeBa.data’ in the current directory.

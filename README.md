# cread2

Software tools for regulatory sequence analysis

CREAD: Comprehensive Regulatory Element Analysis and Discovery.
The CREAD software package is a computational pipeline for understanding
how DNA and RNA elements participate in regulating gene expression.
CREAD uses Pattern-Feature-Model framework machine learning and pattern 
visualization to to identify motifs involved in transcriptional
regulation. CREAD also includes code libraries to facilitate
implementation of new tools. 

Building and Installing
=======================

You may download the latest version of CREAD from: http://smithlabresearch.org/software/cread
This software has been designed to run in a UNIX-like environment.

* Step 0

  This software package requires a functioning installation of the GNU
  Scientific Library (GSL). If you don't already have this installed, you
  will need to download and install it from http://www.gnu.org/software/gsl/

  If gsl is not installed in the default path,
  ```
  export CPATH=/path_to_my_gsl/include
  export LIBRARY_PATH=/path_to_my_gsl/lib
  ```
  will add search paths for compiling and linking.

* Step 1
  
  To build the binaries type the following where '>' is your prompt and the CWD is the root of the distribution:
	
	> make

  This will create all of the executables needed for the CREAD pipleine.

Contacts and Bug Reports
========================

Andrew D. Smith
andrewds@usc.edu

Copyright and License Information
=================================

Copyright (C) 2005-2016
University of Southern California,
Andrew D. Smith

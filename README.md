# cread

Software tools for regulatory sequence analysis

CREAD: Comprehensive Regulatory Element Analysis and Discovery. The
CREAD software package is a computational pipeline for understanding
how DNA and RNA elements participate in regulating gene expression.
CREAD uses Pattern-Feature-Model framework for setting up analyses
with machine learning and pattern visualization to identify motifs
involved in transcriptional regulation. CREAD also includes code
libraries to facilitate implementation of new tools.

Building and Installing
=======================

This software has been designed to run in a UNIX-like environment.

* Step 0

  This software package requires a functioning installation of the GNU
  Scientific Library (GSL). If you don't already have this installed,
  you can download it [here](http://www.gnu.org/software/gsl/). You can
  also get it from apt, conda or brew.

  If GSL is not installed in the default path,
  ```
  export CPATH=/path_to_my_gsl/include
  export LIBRARY_PATH=/path_to_my_gsl/lib
  ```
  will add search paths for compiling and linking.

* Step 1

  To build the binaries type the following from the root of the
  source tree:
  ```console
  make
  make install
  ```
  This will create all of the binaries and then create a `bin`
  directory, in the current directory, containing all the binaries.

Contacts and Bug Reports
========================

Andrew D. Smith
andrewds@usc.edu

Copyright and License Information
=================================

Copyright (C) 2005-2023
Andrew D. Smith

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

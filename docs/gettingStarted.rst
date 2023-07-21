===============
Getting started
===============

SNEEP is fast method to identify regulatory non-coding SNPs (rSNPs) that modify the binding sites of Transcription Factors (TFs) for large collections of SNPs provided by the user. SNEEP combines our statistical approach introduced in our preprint `A statistical approach to identify regulatory DNA variants <https://www.biorxiv.org/content/10.1101/2023.01.31.526404v1>`_ with epigenetic data and additional functionalities.

A graphical summary of SNEEP is shown below:

TODO: add figure

Installation 
==============
We provide a bioconda package to install the main functionality of our approach. Therefore an installation of  Bioconda `here <https://bioconda.github.io/>`_ is requiered. 

TODO: add command

If you want to directly install SNEEP, please clone our `GitHub repository <https://github.com/SchulzLab/SNEEP/>`_ and make sure that the following software is available on your machine: 

- c++11 
- g++ (9.3.0)
- python3.x
- bedtools (v2.27.1)
- openmp

To build SNEEP, run the following commands: 

.. code-block:: console

  cd SNEEP/src/
  make


Please add the path to our software (SNEEP/src) to our PATH environment (otherwise internally called scrips might not be found)

We tested the code and the Makefile only on a linux machine. 

Testing your installation 
==========================

We provide a test script to verify if your installation worked. To download the test data and scripts, please clone the latest version of our GitHub repository:

.. code-block:: console

  clone git@github.com:SchulzLab/SNEEP.git

and our `Zenodo repository <https://doi.org/10.5281/zenodo.4892591>`_. Unzip the files. 

Additional a reference genome in fasta format is requiered. The different chromosomes within the file must be named as chr1, chr2. An example is shown below:


.. example::
  >chr1
  ATCGGGTCA…
  >chr2
  TTTGAGACCAT…


To run our tests, please direct into the SNEEP folder downloaded from GitHub and perform 

.. code-block:: console

  bash runTests.sh <pathToGenome> <pathToZenodoDir>







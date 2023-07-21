===============
Getting started
===============

SNEEP is fast method to identify regulatory non-coding SNPs (rSNPs) that modify the binding sites of Transcription Factors (TFs) for large collections of SNPs provided by the user. SNEEP combines our statistical approach introduced in our preprint is `A statistical approach to identify regulatory DNA variants <https://www.biorxiv.org/content/10.1101/2023.01.31.526404v1>`_.

A summary of the functionalities is shown below:

TODO: add figure

Installation of SNEEP
=======================
We provide a bioconda package to install the main functionality of our approach. Therefore an installation of  Bioconda `here <https://bioconda.github.io/>`_ is requiered. 

TODO: add command

If you want to directly install SNEEP, please clone our `GitHub repository <https://github.com/SchulzLab/SNEEP/>`_ and make sure that the following software is available on your machine: 

- c++11 
- g++ (9.3.0)
- python3.x
- bedtools (v2.27.1)
- openmp

To build SNEEP, run the following commands: 


  cd SNEEP/src/
  make


We tested the code and the Makefile only on a linux machine. 

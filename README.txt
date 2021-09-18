PhyloSignare_v0.1.1  
(Copyright 2021, Authors and Temple University; see license below)

Updated March 28, 2021
==================

The PhyloSignare identifies mutational signatures for each branch of a clone phylogeny [1]. PhyloSignare is written in Python and work smoothly on Windows at present.  It is under continuous development, so anticipate releasing a version for Linux in the near future.  You are free to download, modify, and expand this program under a permissive license similar to the BSD 2-Clause License (see below). PhyloSignare has a few standard dependencies that are noted below.


Dependencies
==================
1. R (version 3.5.2 was tested)
 Please make sure Rscript command is functional.
 R dependencies (signature refitting methods): SignatureEstimation (QP), deconstructSigs (dSig), or MutationalPatterns (MutPat) 
 Note: These dependencies can be obtained from the links listed below.
    SignatureEstimate: https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi\#signatureestimation
    deconstructSigs: https://github.com/raerose01/deconstructSigs
    MutationalPatterns: https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html

2. python (version 3)
 python dependencies (optional):
    pydot
    graphviz
    Note: If the installation of these python packages is not easy, you may want to use Anaconda (https://www.anaconda.com). 


How to prepare input files
==================
1. Map all mutations on a clone phylogeny and make mutation count tables for each branch of a clone phylogeny.
Each mutation needs to be classified into the conventional 96 mutation types (https://cancer.sanger.ac.uk/cosmic/signatures/SBS/). A single file should be created for each branch, and an observed count of each mutation type needs to be listed by using the following format: 
  A[C>A]A,1
  A[C>A]C,4
  A[C>A]G,1
  A[C>A]T,0
  A[C>G]A,2
  A[C>G]C,2
  ...
Example file: Example\Test-A.csv

2. List the mutation count tables and the topology of clone phylogeny.
This file needs to contain the topology of a clone phylogeny and information of branch ID and its corresponding mutation count table, which is prepared above (step 1). To describe the topology of a clone phylogeny, all ancestor and its direct descendant branch pairs should be listed by separating with “->.” For example, B1->B2 means B1 is ancestral branch, and B2 is its direct descendant branch. 
For example,
  #BranchID	File 
  B4	A.csv
  B5	B.csv
  …
  #Tree
  B1->B2
  B5->B3
  ..
Example file: Example\Test.input

3. Make a control file.
A control file should contain the information of a signature refitting method to be used (QP, dSig, MutPat, or MutCon; see the section of R dependencies for the detail), the version of COSMIC signatures (2 or 3; https://cancer.sanger.ac.uk/cosmic), list of expected signatures in your dataset, and ID for the output file name. 
For example, 
  BaseMethod	QP	#QP, dSig, or MutPat
  Signature	2	#2 or 3. 2 for V2 COSMIC signatures and 3 for V3.
  Signature List	all	#all for using all signatures or list of signatures, e.g., S1,S2,S3,S12,S14 for COSMIC V2 and SBS1,SBS2,SBS3 for COSMIC V3
  Signature ID	Test
	
Example file: Control.txt

4. (optional) Edit a color code file.
Color.txt has the color information for each signature. Please assign colors for signatures you like to visualize in the output figure. If you use this option, please install the python dependencies (pydot and graphviz).


How to run PhyloSignare
==================
Open command prompt and run the command, 
python phylosignare.py [your input file] [your control file] 
or
python phylosignare.py [your input file] [your control file] Color.txt

For example: 
To perform the example data analysis, try:
python phylosignare.py Example\Test.input Control.txt 

If you like to generate a figure, in which all signatures are mapped on a clone phylogeny, please use the following command:
python phylosignare.py Example\Test.input Control.txt Color.txt
If you use this option, please install the python dependencies (pydot and graphviz).
For this example dataset, the computation time will be a few minutes. 

Output file
==================
A new directory will be produced in the folder that contains the input file, in which all the output files are found there. 

1. PhyloSignare inference 
In PhyloSignare.txt, signatures that are detected are listed for each branch.

2. PhyloSignare inference (figure)
If color code is provided (e.g., Color.txt) is provided, a figure, in which all signatures are mapped on a clone phylogeny will be produced.

3. Summary file
Mutation count files used and tree topology given are listed in Summary.txt. 

4. Mutation count table
Mutation count tables that are provided are saved, e.g., table for branch ID B5 will be named as B5_MutCount.csv. 

 
Datasets
==================
All datasets used in ref [1] are found at the directory of input_files.

Reference:
[1] Sayaka Miura, Tracy Vu, Jiyeong Choi, Jeffrey P. Townsend, and Sudhir Kumar, Mutational processes in somatic cancer cell populations (2021) Under Review

--------
Copyright 2021, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

--------------------------------
Downloaded from www.pogo-fea.com
--------------------------------

All content Copyright (c) 2013 Peter Huthwaite.

Directories
------------
blocker - the aligned partitioner
blockgreedy - the greedy partitioner
common - some common code, including the parameters (e.g. block sizes)
inputGenerate - code to convert Abaqus input files to Pogo input files
solver - the pogo solver itself

Note:
The contents of parser/ and common/eigen3/ are open source 
components from other authors, and I do not claim to have 
written this code. The parser is used to process text in 
the Abaqus input file, and eigen is a matrix multiplication 
library I use during pre-processing to calculate element 
stiffness.

To compile
------------
The makefiles should 'just work'. None of the compilation is 
particularly complex, so it should be straightforward to 
identify what is going wrong if necessary. Note that depending on 
your setup, you may have to change options at the top of the 
makefiles (e.g. the gcc version etc. so it plays nicely with CUDA).

To use
---------
From an Abaqus input file 'test.inp':

$./abqConvert test.inp

This should produce a file test.pogo-inp. Then run:

$./blocker test

(or blockgreedy). This should take test.pogo-inp and generate 
test.pogo-block. Then the solver can be run:

$./pogo test

(or pogo64) which should take both test.pogo-inp and test.pogo-block. 



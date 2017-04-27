This is a MATLAB-implementation of the asymmetric version of gplmDCA.
The program uses 'minFunc' by Mark Schmidt (contained in the folder '3rd_party_code'), which can be found at http://www.di.ens.fr/~mschmidt/Software/minFunc.html.
 
Inputs to gplmDCA: 
	fastafile:
	Alignment file in FASTA format. Inserts, which should be represented by '.' and lower-case letters (as is standard in the Pfam download), are removed automatically by the program.

	outputfile:
	This becomes a file with N(N-1)/2 rows (N=domain length), each with three entries: residue i, residue j, and interaction score for the pair (i,j).

	lambda_h:
	Field-regularization strength (typical value: 0.01).

	lambda_J:
	Coupling-regularization strength (typical value: 0.01 - this is tested for N in the range 50-100 or so, and the optimal value may be different for longer domains).

	lambda_chi:
	Gap-regularization strength (typical value 0.001).

	reweighting_threshold:
	Required fraction of nonidentical AA for two sequences to be counted as independent (typical values: 0.1-0.3).
	Note that this is not the threshold 'x' as defined in the paper, but 1-x.
 	
	nr_of_cores:
	The number of processors to use on the local machine. 
	If this argument is >1, the program calls functions from MATLAB's Parallel Computing Toolbox.

	M:
	Maximal gap length for gap parameters. 
	If set to -1, it will be chosen automatically, otherwise only gaps with length [1, M] will be affected. Unless
	you have a good reason to do so, we recommend setting M to -1.
	

Typical call (on a quad-core machine):
	plmDCA_asymmetric('alignment.fas','alignment.gplmdca',0.01,0.01,0.1,4)


---------------------------------------------
Copyright conditions for gplmDCA:
	Copyright 2014 - Christoph Feinauer and Marcin J. Skwark (christophfeinauer@gmail.com, marcin@skwark.pl)
	Copyright 2012 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
	
	All rights reserved

	Permission is granted for anyone to copy, use, or modify this
	software for any uncommercial purposes, provided this copyright 
	notice is retained, and note is made of any changes that have 
	been made. This software is distributed without any warranty, 
	express or implied. In no event shall the author or contributors be 
	liable for any damage arising out of the use of this software.

	The publication of research using this software, modified or not, must include an 
	appropriate citation to:
	
	C. Feinauer, M.J. Skwark, A. Pagnani, E. Aurell, Improving Contact Prediction 
	along Three Dimensions. PLoS Comput Biol 10(10): e1003847. 
	doi: 10.1371/journal.pcbi.1003847

	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013) 
	doi: 10.1103/PhysRevE.87.012707

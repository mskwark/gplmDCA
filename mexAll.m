fprintf('Compiling g_rC.c\n');
mex -outdir functions functions/g_rC.c

fprintf('Compiling calc_inverse_weights.c\n');
mex -outdir functions functions/calc_inverse_weights.c

fprintf('Compiling gapCount.c\n');
mex -outdir functions functions/gapCount.c

fprintf('Compiling gapMat.c\n');
mex -outdir functions functions/gapMat.c

fprintf('Compiling lbfgsAddC.c in minFunc\n');
mex -outdir 3rd_party_code/minFunc/ 3rd_party_code/minFunc/lbfgsAddC.c

fprintf('Compiling lbfgsC.c in minFunc\n');
mex -outdir 3rd_party_code/minFunc/ 3rd_party_code/minFunc/lbfgsC.c

fprintf('Compiling lbfgsProdC.c in minFunc\n');
mex -outdir 3rd_party_code/minFunc/ 3rd_party_code/minFunc/lbfgsProdC.c

fprintf('Compiling mcholC.c in minFunc\n');
mex -outdir 3rd_party_code/minFunc/ 3rd_party_code/minFunc/mcholC.c

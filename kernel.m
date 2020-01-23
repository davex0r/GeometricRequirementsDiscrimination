function p = kernel(L)
% Compute the normalized (to sum 1) basis for the kernel of L, does not
% check if L is a valid Laplacian matrix, use with care!

 p = null(L)/sum(null(L));
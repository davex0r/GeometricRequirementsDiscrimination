function L = makeLaplacian(L)
% Impose the laplacian condiiton on a matrix, does not check if the matrix
% is square or if all the netries are positive, use with care. 

L  = L.*(ones(size(L))-eye(size(L)));
L = L+diag(-sum(L,1));


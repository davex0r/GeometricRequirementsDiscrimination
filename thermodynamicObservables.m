function [S_i, S_e] = thermodynamicObservables(L,p)
% Compute the entropy production Si and the entropy flow Se following
% Schankenberg(1974) and Rao,Peliti(2015).

J = L'.*p;
tmp = (triu(J,1)-tril(J,-1)').*log(triu(J,1)./tril(J,-1)');
tmp2 = L';
tmp1 = (triu(J,1)-tril(J,-1)').*log(triu(tmp2,1)./tril(tmp2,-1)');
S_i = 1/2*nansum(tmp(:));
S_e = -1/2*nansum(tmp1(:));
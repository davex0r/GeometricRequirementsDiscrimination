function [error_rate, L, p, G]= fixedKineticLadder(s,vars,gamma,delta)
% Irreversible kinetic ladder graph following Murugan, Huse and Leibler
% (2012).  Error rate, as well as teh Laplacian Matric (L) and the steady
% state distribution (p) are returned.  "s" is the number of nodes and must
% be 9+4N (9 is 1 loop)

f = vars(1);
b = vars(2);
u = vars(3);
d = vars(4);

L = zeros(s);

L(1,2:3) = [b b];
L(2,[4 6]) = [d*exp(delta) b];
L(3,[5 7]) = [d*exp(gamma) b];
L(4,[1 2]) = [f u*exp(delta)];
L(5,[1 3]) = [f u];

n = (s-5-4)/4;
for i = 1:n
    ii = 6+(i-1)*4;
    L(ii,[ii+2 ii+4]) = [d*exp(delta) b];
    L(ii+1,[ii+3 ii+5]) = [d*exp(gamma) b];
    L(ii+2,[ii ii-2]) = [u*exp(delta) f];
    L(ii+3,[ii+1 ii-1]) = [u f];
end

L(s,[s-2 s-4]) = [u f];
L(s-1,[s-3 s-5]) = [u*exp(delta) f];
L(s-2,s) = d*exp(gamma);
L(s-3,s-1) = d*exp(delta);


    
% L = L+diag(-sum(L),0);
L = makeLaplacian(L);
G = digraph(L','omitselfloops');
p = kernel(L);
if all(p>=0)
    error_rate = p(end)/p(end-1);
else
    error_rate = NaN;
end

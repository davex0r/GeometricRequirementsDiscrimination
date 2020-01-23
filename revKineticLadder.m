function [error_rate, L, p, G]= revKineticLadder(s,vars,w,epsilon,epsilon_b,epsilon_f,epsilon_u,gamma,delta,F)
% Reversible Ladder model as in Rao and Peliti, J. Stat Mech 2015.
% Error rate, as well as teh Laplacian Matric (L) and the steady
% state distribution (p) are returned.  "s" is the number of nodes and must
% be 9+4N (9 is 1 loop)

w_f = vars(1);
w_b = vars(2);
w_u = vars(3);

krp = w*exp(epsilon);
krm = w;

fp = w_f;
fm = w_f*exp(-epsilon_f);
bp = w_b;
bm = w_b*exp(-epsilon_b);
dr = w_u*exp(delta);
dw = w_u*exp(gamma);
ur = w_u*exp(epsilon_u+delta);
uw = w_u*exp(epsilon_u);

L = zeros(s);

L(1,2:3) = [krm krm];
L(2,[4 6]) = [dr bp];
L(3,[5 7]) = [dw bp];
L(6,2) = bm;
L(7,3) = bm;
L(4,2) = ur;
L(5,3) = uw;
L([2 3],1) = krp;

n = (s-5-4)/4;
for i = 1:n
    ii = 6+(i-1)*4;
    L(ii,[ii+2 ii+4]) = [dr bp];
    L(ii+4,ii) = bm;
    L(ii+1,[ii+3 ii+5]) = [dw bp];
    L(ii+5,ii+1) = bm;
    L(ii+2,[ii ii-2]) = [ur fp];
    L(ii-2,ii+2) = fm;
    L(ii+3,[ii+1 ii-1]) = [uw fp];
    L(ii-1,ii+3) = fm;
end

L(s,[s-2 s-4]) = [uw fp];
L(s-4,s) = fm;
L(s-1,[s-3 s-5]) = [ur fp];
L(s-5,s-1) = fm;
L(s-2,s) = dw;
L(s-3,s-1) = dr;

L(1,end-1:end) = L(1,end-1:end)+F;


    
L = L+diag(-sum(L),0);
G = digraph(L','omitselfloops');
p = kernel(L);
if all(p>=0)
    error_rate = p(end)/(p(end-1));
else
    error_rate = NaN;
end

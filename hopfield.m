function [error_rate, L, p, G] = hopfield(vars,w,epsilon,gamma,delta,delta_p,F)
% Generate the Laplacian Matrix (L), and compute the error rate and steady
% state distribution (p) for the Hopfield Ninio model with Kramer's form
% rate contants

w_i = (vars(1));
w_p = (vars(2));
epsilon_i = (vars(3));
epsilon_p = (vars(4));

krp = w*exp(delta+epsilon);
krm = w*exp(delta);
kwp = w*exp(epsilon);
kwm = w*exp(gamma);

hrp = w_i*exp(epsilon_i);
hwp = w_i*exp(epsilon_i);

hrm = w_i;
hwm = w_i;

Krp = w_p*exp(-delta_p);
Kwp = w_p;

Krm = w_p*exp(epsilon_p-delta_p);
Kwm = w_p*exp(epsilon_p+gamma);


L = [...
    [-(krp+Krp+kwp+Kwp) krm Krm+F kwm Kwm+F];...
    [krp -(krm+hrp) hrm 0 0];...
    [Krp hrp -(Krm+F+hrm) 0 0];...
    [kwp 0 0 -(kwm+hwp) hwm];...
    [Kwp 0 0 hwp -(Kwm+F+hwm)]];

p = kernel(L);
if all(p>=0)
    error_rate = p(5)/(p(3));
else
    error_rate = 1;
end

G = digraph(L','omitselfloops');
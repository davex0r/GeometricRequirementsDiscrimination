function [error_rate, L, p, G] = double_hopfield(vars,w,epsilon,gamma,delta,delta_p,F)
% Special Hopfield Ninio graph for product tuning. In this graph, delta favors one
% product while gamma favors the other.  It is otherwise identical to
% hopfield.m

w_i = (vars(1));
w_p = (vars(2));
epsilon_i = (vars(3));
epsilon_p = (vars(4));

krp = w*exp(epsilon);
krm = w;
kwp = w*exp(epsilon+delta);
kwm = w*exp(gamma+delta);

hrp = w_i*exp(epsilon_i);
hwp = w_i*exp(epsilon_i);

hrm = w_i;
hwm = w_i;

Krp = w_p;
Kwp = w_p*exp(-delta_p);

Krm = w_p*exp(epsilon_p);
Kwm = w_p*exp(epsilon_p+gamma-delta_p);


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
%% Numeric Simulation (Hopfield,Gamma)
clear; clc
clf

w = 1;
epsilon = 10;
gamma = 1;

rng(1)
iMax = 1000;

vars = [exp(0.5*randn(iMax,2)) 2*(randn(iMax,2))];
F = 0;
delta = 0;
delta_p = 0;
L_out = nan(5,5,iMax);
error_rate = nan(1,iMax);
S_i = nan(1,iMax); S_e = nan(1,iMax); non_orth = nan(1,iMax);

for i = 1:iMax
    [error_rate(i), L, p, ~] = hopfield(vars(i,:),w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-6) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        dno = norm(S'*S-eye(size(S,2)),'fro');
        non_orth(i) = (1-dno);
        cyclic_balance(i,1) = (L(1,2)*L(2,3)*L(3,1))/(L(1,3)*L(3,2)*L(2,1));
    end
end
scatter(non_orth,log(error_rate),10,log10(S_i))
xlims = get(gca,'xlim');
set(gca,'fontsize',14,'fontweight','bold')
xlabel('\Theta')
ylabel('ln(\xi)')
line(xlims,([(-(gamma)) (-(gamma))]),'color','g')
line(xlims,([(-(2*gamma)) (-(2*gamma))]),'color','r')
%% Numeric Simulation (Hopfield,Delta)
clear; clc
clf

w = 1;
epsilon = 10;
gamma = 0;

rng(1)
iMax = 1000;

vars = [exp(0.5*randn(iMax,2)) 2*(randn(iMax,2))];
F = 0;
delta = 1;
delta_p = 1;
L_out = nan(5,5,iMax);
error_rate = nan(1,iMax);
S_i = nan(1,iMax); S_e = nan(1,iMax); non_orth = nan(1,iMax);

for i = 1:iMax
    [error_rate(i), L, p, ~] = hopfield(vars(i,:),w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-6) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        dno = norm(S'*S-eye(size(S,2)),'fro');
        non_orth(i) = (1-dno);
    end
end
scatter(non_orth,log(error_rate),10,log10(S_i))
xlims = get(gca,'xlim');
set(gca,'fontsize',14,'fontweight','bold')
xlabel('\Theta')
ylabel('ln(\xi)')
line(xlims,([(-(delta)) (-(delta))]),'color','g')
line(xlims,([(-(2*delta)) (-(2*delta))]),'color','r')
%% Orthogonality and Dissipation Minimize and Maximize (Hopfield, Gamma)
clear
clc
tmp = [0.55    0.7318   -3.3256    4.2245];
w = 1;
epsilon = 10;
gamma = 1;
delta = 0;
delta_p = 0;
F = 0;

xlims = [-10 10];
mp = linspace(xlims(1),xlims(2),100);

for i = 1:length(mp)
    tmp1 = tmp;
    tmp1(3) = mp(i);
    [error_rate(i), L, p, ~] = hopfield(tmp1,w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-4) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        dno = norm(S'*S-eye(size(S,2)),'fro');
        non_orth(i) = (1-dno);
    end
end

[old_co, new_co] = init_color;
yyaxis left
plot(mp*tmp(1),(non_orth),'linewidth',2,'color',old_co(2,:))
set(gca,'ycolor',old_co(2,:))
ylabel('\Theta')
set(gca,'ylim',[-1 1],'ytick',[-1 0 1])
hold on
yyaxis right
plot(mp*tmp(1),((log(error_rate))),'-','linewidth',2,'color',old_co(1,:));
set(gca,'ycolor',old_co(1,:))
[~, indx] = min(error_rate);

xlabel('log(m^\prime)')
diss = log10(S_i);
diss = ((diss-min(diss))/range(log10(S_i))*2)-2;
set(gca,'ytick',[-2 -1])
ylabel('log(\xi)')

p1 = plot(mp*tmp(1),diss,'-','linewidth',2,'color',old_co(4,:));
ylims = get(gca,'ylim');
line([mp(indx)*tmp(1) mp(indx)*tmp(1)],ylims,'color','r')
legend(p1,'Dissipation (scaled)','Location','southeast')
legend('boxoff')
set(gca,'ytick',[-2 -1 0],'ylim',[-2 0],'xlim',xlims*tmp(1),'fontsize',14,'fontweight','bold')

%% Orthogonality and Error (Hopfield, Gamma, analytic S'S elements)
clear
clc
tmp = [0.55    0.7318   -3.3256    4.2245];

w = 1;
epsilon = 10;
gamma = 1;
delta = 0;
delta_p = 0;
F = 0;
[old_co, new_co] = init_color;

iMax = 50;
dno = nan(3,3,iMax);
xlims = [-4 4];
ei = linspace(xlims(1),xlims(2),iMax);
for i = 1:length(ei)
    tmp1 = tmp;
    tmp1(3) = ei(i);
    mu = tmp1(1)*exp(tmp1(3));
    wp = tmp1(2);
    [error_rate(i), L, p, ~] = hopfield(tmp1,w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-4) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
%         dno = norm(S'*S-eye(size(S,2)),'fro');
%         non_orth(i) = (max_orth-dno)/5;
        dno(:,:,i) = abs(S'*S);
        orth(i) = 1-(norm(eye(size(S,2))-dno(:,:,i),'fro'));
        a(i) = sqrt((w^4*exp(2*gamma))/(4*(w^2+w*mu+mu^2)*(w^2*exp(2*gamma)+w*exp(gamma)*mu+mu^2)));
        b0 = (wp*(2*w-mu)+exp(epsilon)*w*(3*w+mu))^2;
        b1 = 4*(3*exp(2*epsilon)*w^2+4*exp(epsilon)*w*wp+3*wp^2)*(w^2+w*mu+mu^2);
        b(i) = sqrt(b0/b1);
        c0 = (3*exp(gamma+epsilon)*w^2+2*exp(gamma)*w*wp+exp(epsilon)*w*mu+wp*mu)^2;
        c1 = (4*(3*exp(2*epsilon)*w^2+4*exp(epsilon)*w*wp+3*wp^2)*(exp(2*gamma)*w^2+exp(gamma)*w*mu+mu^2));
        c(i) = sqrt(c0/c1);
    end
end

var1 = ei;
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on
set(gca,'ycolor',old_co(2,:))
pp(1) = plot((var1),a,'-','color',new_co(1,:));
pp(2) = plot((var1),b,'-','color',new_co(2,:));
pp(3) = plot((var1),c,'-','color',new_co(3,:));
set(gca,'ylim',[-1 1],'ytick',[-1 0 1])
ylabel('\Theta')
yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
xlabel('log(\epsilon_i)')
ylabel(' log(\xi)')
set(gca,'ytick',[-2 -1 0],'ylim',[-2 0],'fontsize',14,'fontweight','bold')
legend(pp,{'<2,3>','<1,2>','<1,3>'},'Location','northwest')
text(-0.3,-0.1,'\gamma=1, \delta=0','fontsize',18)
%% Orthogonality and Error (Hopfield, Delta, analytic S'S elements)
clear
clc
tmp = [2.2699    0.8982    6.2196    2.5553];
[old_co, new_co] = init_color;

iMax = 50;
dno = nan(3,3,iMax);
xlims = [-4 4];
ei = linspace(xlims(1),xlims(2),iMax);
w = 1;
epsilon = 10;
gamma = 0;
delta = 1;
delta_p = 1;
F = 0;
ei = linspace(-4,4,50);
for i = 1:length(ei)
    tmp1 = tmp;
    tmp1(3) = ei(i);
    mu = tmp1(1)*exp(tmp1(3));
    wp = tmp1(2);
    [error_rate(i), L, p, ~] = hopfield(tmp1,w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-4) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        dno(:,:,i) = abs(S'*S);
        orth(i) = 1-(norm(eye(size(S,2))-dno(:,:,i),'fro'));
        a(i) = sqrt((w^4*exp(2*delta))/(4*(w^2+w*mu+mu^2)*(w^2*exp(2*delta)+w*exp(delta)*mu+mu^2)));
        
        b0 = (exp(delta)*(1+exp(delta_p))*w*wp-wp*mu+exp(delta+delta_p+epsilon)*w*(w+2*exp(delta)*w+mu))^2;
        b1 = 4*(exp(2*delta_p+2*epsilon)*w^2+exp(2*(delta+delta_p+epsilon))*w^2+...
            exp(delta_p+epsilon)*w*wp+exp(2*delta_p+epsilon)*w*wp+wp^2+...
            exp(delta_p)*wp^2+exp(2*delta_p)*wp^2+exp(delta+delta_p+epsilon)*w*(wp+exp(delta_p)*(exp(epsilon)*w+wp)))...
            *(exp(2*delta)*w^2+exp(delta)*w*mu+mu^2);
        b(i) = sqrt(b0/b1);
        c0 = (w*(-exp(epsilon)*w-exp(delta+epsilon)*w-wp-exp(-delta_p)*wp)+...
            (exp(epsilon)*w*(-w-mu)+wp*mu))^2;
        c1 = (exp(2*epsilon)*w^2+exp(2*delta+2*epsilon)*w^2+wp^2+exp(-2*delta_p)*wp^2+...
            (exp(epsilon)*w+exp(delta+epsilon)*w+wp+exp(-delta_p)*wp)^2)*...
            (w^2+(w+mu)^2+mu^2);
        c(i) = sqrt(c0/c1);
    end
end

var1 = ei;
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on
set(gca,'ycolor',old_co(2,:))
pp(1) = plot((var1),a,'-','color',new_co(1,:));
pp(2) = plot((var1),b,'-','color',new_co(2,:));
pp(3) = plot((var1),c,'-','color',new_co(3,:));
set(gca,'ylim',[-1 1],'ytick',[-1 0 1])
ylabel('\Theta')
yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(\xi)')
xlabel('log(\epsilon_i)')
text(-0.3,-0.1,'\gamma=0, \delta=1','fontsize',18)
set(gca,'ytick',[-2 -1 0],'ylim',[-2 0],'fontsize',14,'fontweight','bold')
%% Product selection via Dissapation across a single reaction node (Hopfield, Custom Gamma Delta)
clear
clc
tmp = [0.6245    0.1917   -2.6672    5.2397];
w = 1;
epsilon = 5;
gamma = 1;
delta = 3;
delta_p = 0;
F = 0;

xlims = [-2 8];
ei = linspace(xlims(1),xlims(2),50);
for i = 1:length(ei)
    tmp1 = tmp;
    tmp1(3) = ei(i);
    mu = tmp1(1)*exp(tmp1(3));
    wp = tmp1(2);
    [error_rate(i), L, p, ~] = double_hopfield(tmp1,w,epsilon,gamma,delta,delta_p,F);
    if all(L*p<1e-4) &&all(p>=0)
        L_out(:,:,i) = L;
        S = L(:,[1 2 4]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        dno(:,:,i) = abs(S'*S);
        orth(i) = (1-norm(eye(size(S,2))-dno(:,:,i),'fro'))/5;
        a(i) = sqrt((w^4*exp(2*delta))/(4*(w^2+w*mu+mu^2)*(w^2*exp(2*delta)+w*exp(delta)*mu+mu^2)));
        b0 = (exp(delta)*(1+exp(delta_p))*w*wp-wp*mu+exp(delta+delta_p+epsilon)*w*(w+2*exp(delta)*w+mu))^2;
        b1 = 4*(exp(2*delta_p+2*epsilon)*w^2+exp(2*(delta+delta_p+epsilon))*w^2+...
            exp(delta_p+epsilon)*w*wp+exp(2*delta_p+epsilon)*w*wp+wp^2+...
            exp(delta_p)*wp^2+exp(2*delta_p)*wp^2+exp(delta+delta_p+epsilon)*w*(wp+exp(delta_p)*(exp(epsilon)*w+wp)))...
            *(exp(2*delta)*w^2+exp(delta)*w*mu+mu^2);
        b(i) = sqrt(b0/b1);
        c0 = (w*(-exp(epsilon)*w-exp(delta+epsilon)*w-wp-exp(-delta_p)*wp)+...
            (exp(epsilon)*w*(-w-mu)+wp*mu))^2;
        c1 = (exp(2*epsilon)*w^2+exp(2*delta+2*epsilon)*w^2+wp^2+exp(-2*delta_p)*wp^2+...
            (exp(epsilon)*w+exp(delta+epsilon)*w+wp+exp(-delta_p)*wp)^2)*...
            (w^2+(w+mu)^2+mu^2);
        c(i) = sqrt(c0/c1);
    end
end
clf

var1 = ei;

[old_co, new_co] = init_color;
yyaxis left
pp(1) = plot((var1),orth,'-','linewidth',2,'color',old_co(2,:));
set(gca,'ycolor',old_co(2,:),'ylim',[-0.2 0.1]);
ylabel('\Theta')
hold on
diss = log10(S_i);
pp(3) = plot(var1,((diss-min(diss))/range(diss))*0.3-0.2,'-','linewidth',2,'color',old_co(4,:));
yyaxis right
pp(2) = plot((var1),-log(error_rate),'linewidth',2,'color',old_co(1,:));
set(gca,'ycolor',old_co(1,:));
ylabel('log(p_\gamma/p_\delta)')
line(xlims,[0 0],'color','r','linestyle','--')
xlabel('(\epsilon_i)')
set(gca,'xlim',xlims,'ytick',[-2 -1 0 1 2])
set(gca,'linewidth',2,'fontsize',14,'fontweight','bold')
legend(pp(3),'Dissipation (scaled)','location','south')
legend('boxoff')
%% Orthogonality and Error (Ladder, Gamma, analytic S'S elements)
clear
clc
[old_co, new_co] = init_color;
iMax = 100;
tmp = [0.1 2 0.1 2];
xlims = [-10 2];
x = logspace(xlims(1),xlims(2),iMax);
s = 17;
for i = 1:iMax
    tmp1 = tmp;
    tmp1(4) = x(i);
    [error_rate(i),L,p,~] = fixedKineticLadder(s,tmp1,1,0);
    S = L(:,1:end-2);
    S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
    dno(:,:,i) = abs(S'*S);
    orth(i) = 1-norm(eye(size(S,2))-dno(:,:,i),'fro');
end
var1 = log10(x);
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on
set(gca,'ycolor',old_co(2,:),'ylim',[-3 -1],'ytick',[-3 -2 -1])
ylabel('\Theta')
yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(\xi)')
xlabel('log(\omega_u)')
set(gca,'ytick',[-4 -3 -2 -1],'ylim',[-4 -1])
text(-0.8,-1.2,'\gamma=1, \delta=0','fontsize',18)
%% Orthogonality and Error (Ladder, Delta, analytic S'S elements)
[old_co, new_co] = init_color;
iMax = 100;
tmp = [2 0.1 3 0.1];
xlims = [-5 2];
x = logspace(xlims(1),xlims(2),iMax);
s = 17;
for i = 1:iMax
    tmp1 = tmp;
    tmp1(4) = x(i);
    [error_rate(i),L,p,~] = fixedKineticLadder(s,tmp1,0,1);
    S = L(:,1:end-2);
    S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
    dno(:,:,i) = abs(S'*S);
    orth(i) = 1-norm(eye(size(S,2))-dno(:,:,i),'fro');
end
var1 = log10(x);
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on
set(gca,'ycolor',old_co(2,:),'ylim',[-3 -1],'ytick',[-3 -2 -1])
ylabel('\Theta')
yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(\xi)')
xlabel('log(\omega_u)')
set(gca,'ytick',[-3 -2 -1 0],'ylim',[-3 0.5])
text(-0.8,0.25,'\gamma=0, \delta=1','fontsize',18)
%% Product selection via Dissapation across a single reaction node (Ladder, Custom Gamma Delta)
clear 
clc
iMax = 50;
[old_co, new_co] = init_color;
tmp = [  0.0468    0.6582    5.2381];
xlims = [-3 2];
x = logspace(xlims(1),xlims(2),iMax);
s = 17;
gamma = 1;
delta = 2;
epsilon = 10;
w = 1;
epsilon_b = 8;
epsilon_f = 7;
epsilon_u = 7;
F = 0;
for i = 1:iMax
    tmp1 = tmp;
    tmp1(3) = x(i);
    [error_rate(i),L,p,~] = doubleThermoLadder(s,tmp1,w,epsilon,epsilon_b,epsilon_f,epsilon_u,gamma,delta,F);
    S = L(:,1:end-2);
    S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
    [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
    dno(:,:,i) = abs(S'*S);
    orth(i) = 1-norm(eye(size(S,2))-dno(:,:,i),'fro');
end
var1 = -log10(x);
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on
diss = log(S_i);
pp(3) = plot((var1),(diss-min(diss))/range(diss)-3,'-','linewidth',2,'color',old_co(4,:));
set(gca,'ycolor',old_co(2,:),'ylim',[-3 -2],'ytick',-3:0.5:-2)
ylabel('\Theta')

yyaxis right
plot((var1),-log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(p_{\gamma}/p_{\delta})')
xlabel('log(\eta)')
set(gca,'linewidth',2,'fontsize',14,'fontweight','bold')
legend(pp(3),'Dissipation','Location','northeast')
legend('boxoff')
line(-fliplr(xlims),[0 0],'color','r','linestyle','--')
legend(pp(3),'Dissipation','Location','southwest')
legend('boxoff')
%% Orthogonality and Dissipation Minimize and Maximize (Ladder, Gamma, different limits)
clear 
clc
iMax = 50;
[old_co, new_co] = init_color;
tmp = [  0.0874    2.3565   15.3333];
xlims = [-5 5];
x = logspace(xlims(1),xlims(2),iMax);
s = 17;

gamma = 1;
delta = 0;
epsilon = 3;
w = 1;
epsilon_b = 3;
epsilon_f = 3;
epsilon_u = 3;
F = 0;

for i = 1:iMax
    tmp1 = tmp;
    tmp1(1) = x(i);
    [error_rate(i),L,p,~] = revKineticLadder(s,tmp1,w,epsilon,epsilon_b,epsilon_f,epsilon_u,gamma,delta,F);
    S = L(:,1:end-2);
    [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
    S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
    dno(:,:,i) = abs(S'*S);
    orth(i) = 1-norm(eye(size(S,2))-dno(:,:,i),'fro');
end
[~, indx] = min(error_rate);
var1 = log10(x);
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on

set(gca,'ycolor',old_co(2,:),'ylim',[-2.75 -1.75],'ytick',[-2.5 -2])
ylabel('\Theta')

yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(\xi)')
xlabel('log(\omega_f)')

diss = log(S_i);
pp(3) = plot((var1),(diss-min(diss))/range(diss)*2.75-4,'-','linewidth',2,'color',old_co(4,:));
set(gca,'ycolor',old_co(1,:),'ylim',[-4 -0.8],'ytick',[-4 -3 -2 -1])
ylims = get(gca,'ylim');
line([var1(indx) var1(indx)],ylims,'color','r')
legend(pp(3),'Dissipation','Location','west')
legend('boxoff')
%% Orthogonality and Dissipation Minimize and Maximize (Ladder, Gamma, same limits)
clear 
clc
iMax = 50;
[old_co, new_co] = init_color;
tmp = [  0.0874    2.3565   15.3333];
xlims = [-5 5];
x = logspace(xlims(1),xlims(2),iMax);
s = 17;

gamma = 1;
delta = 0;
epsilon = 3;
w = 1;
epsilon_b = 3;
epsilon_f = 3;
epsilon_u = 3;
F = 0;

for i = 1:iMax
    tmp1 = tmp;
    tmp1(3) = x(i);
    [error_rate(i),L,p,~] = revKineticLadder(s,tmp1,w,epsilon,epsilon_b,epsilon_f,epsilon_u,gamma,delta,F);
    S = L(:,1:end-2);
    [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
    S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
    dno(:,:,i) = abs(S'*S);
    orth(i) = 1-norm(eye(size(S,2))-dno(:,:,i),'fro');
end
[~, indx] = min(error_rate);
var1 = log10(x);
yyaxis left
plot((var1),orth,'linewidth',2,'color',old_co(2,:))
hold on

set(gca,'ycolor',old_co(2,:),'ylim',[-2.75 -1.75],'ytick',10.5:0.5:12)
ylabel('\Theta')

yyaxis right
plot((var1),log(error_rate),'linewidth',2,'color',old_co(1,:))
set(gca,'ycolor',old_co(1,:))
ylabel(' log(\xi)')
xlabel('log(\omega_u)')

ylims = get(gca,'ylim');
diss = log(S_i);
pp(3) = plot((var1),(diss-min(diss))/range(diss)*2.75-4,'-','linewidth',2,'color',old_co(4,:));
set(gca,'ycolor',old_co(1,:),'ylim',[-4 -1],'ytick',[-3 -2 -1])

line([var1(indx) var1(indx)],ylims,'color','r')
legend(pp(3),'Dissipation','Location','east')
legend('boxoff')
%% Orthogonality Increases as we move toward fully connected graph
clear; clc;
n = 16;
combs = nchoosek(1:n,2);
sparsities = 0.1:0.1:0.9;
nBreaks = floor((1-sparsities)*size(combs,1));
for m = 1:length(sparsities)+2
    for i = 1:1000
        L = exp(randn(n)/3);
%         energies = randn(n,1)/3;
%         L = exp(-(energies-energies'));
        if m==2
            tmp_s = size(L,1);
            tmp_o = sqrt(tmp_s);
            A = diag(ones(tmp_s-tmp_o,1),tmp_o)+diag(ones(tmp_s-tmp_o,1),-tmp_o);
            tmp1 = diag(A,1);
            for k = 1:length(tmp1)
                if rem(k,sqrt(tmp_s))~=0
                    tmp1(k) = 1;
                end
            end
            A = A+diag(tmp1,1)+diag(tmp1,-1);
            
            L = L.*A;
        end
        
        if m>=3
            indx = randi(size(combs,1),nBreaks(m-2));
            for j = 1:length(indx)
                L(combs(indx(j),1),combs(indx(j),2)) = 0;
                L(combs(indx(j),2),combs(indx(j),1)) = 0;
            end
        end
        L = makeLaplacian(L);
        S = L(:,[2:end-1]);
        S = S./(repmat(sum(S.^2),size(S,1),1).^(1/2));
        p = kernel(L);
        [S_i(i), S_e(i)] = thermodynamicObservables(L,p);
        
        
        gram_matrix = S'*S;
        off_diag_elems(:,i) = gram_matrix(triu(true(n-2),1));
        max_orth = norm(eye(n-2)-ones(n-2),'fro');
        omega(i) = 1-norm(eye(n-2)-S'*S,'fro');
        omega_run(m) = nanmean(omega/n);
        omega_run_sd(m) = nanstd(omega/n);
    end
    figure(1)
histogram(omega/n,'Normalization','pdf'); xlabel('Orthogonality per Node (\Omega/n)')
hold on
legend({'Fully Connected','Grid','Sparse'})
end
figure
errorbar([1 0 sparsities],omega_run,omega_run_sd,'s','linewidth',2)
xlabel('Fraction of Connections (0=Grid)')
ylabel('\Theta')
set(gca,'fontsize',14,'fontweight','bold')

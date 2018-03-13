
clear;clc;close all;
%perfectly sparse signals
%seed = round(rand()*100);
%disp(['seed=', num2str(seed)])

rand(1)

m = 500;
Na = 500;
Nb = 500;
N = Na + Nb;
k = 700; %sparsity

%sensing matrix
%A = normc(normrnd(1, 1, [m, Na]));
%B = normc(normrnd(1, 1, [m, Nb]));
%A = normc(rand(m, Na));
%B = normc(rand(m, Nb));

Atemp = dctmtx(m);
Atemp = Atemp(1:m, 1:Na);
Btemp = eye(m, Nb);
Dtemp = [Atemp Btemp];
Dtemp = Dtemp(:,randperm(Na+Nb));

A = normc(Dtemp(:,1:Na));
B = normc(Dtemp(:,Na+1:Na+Nb));

D = [ A B ];

%% Values- Original %%
mu_a = coherence(A);
mu_b = coherence(B);
mu_larger  = max(mu_a, mu_b);
mu_smaller = min(mu_a, mu_b);
    
mu_m = mutualcoherence(A,B);
mu_d = coherence(D);

mu_d_alt = (0.5 * (mu_larger*(k-2) + k*sqrt(mu_larger^2 + mu_m^2)))/(k-1);
mu_d_eff = min(mu_d, mu_d_alt);

k_bound1 = 2*(1 + mu_larger)/(mu_larger + 2*mu_d + sqrt(mu_larger^2 + mu_m^2));
k_bound2 = (1+mu_d)/(2*mu_d);
k_bound = max(k_bound1, k_bound2);

deltahat1 = 0.5 * (mu_larger*(k-2) + k*sqrt(mu_larger^2 + mu_m^2));
deltahat2 = mu_d * (k-1);
deltahat = min(deltahat1, deltahat2);
%deltahat =0;

paramlist = {'Na', 'Nb', 'mu_a', 'mu_b', 'mu_m', 'mu_d', 'mu_d_alt', 'mu_d_eff', 'k_bound', 'deltahat'};
valueset_original = [Na, Nb, mu_a, mu_b, mu_m, mu_d, mu_d_alt, mu_d_eff, k_bound, deltahat];

%% Optimally rearranged & split  (2-way)
tic
G = abs(D' * D); %coherence matrix
G(logical(eye(size(G)))) = 0;

setA = [];
setB = [];
unassigned = 1:N;

[~,I] = max(G(:));
[u, v] = ind2sub(size(G),I);

setA = [setA; u];
setB = [setB; v];
assigned = [setA; setB];
unassigned = setdiff(unassigned, assigned);
G(u,v) = 0;
G(v,u) = 0;

numassigned = 2;
while max(G(:)) > 0
    u=1; v=1;
    maxwt = 0;
    uindex = 0;
    
    for i = 1:length(assigned)
        ku = assigned(i);
        [kmaxwt, kv] = max(G(ku,:));
        
        if kmaxwt > maxwt
            u = ku;
            v = kv;
            maxwt = kmaxwt;
            uindex = i;
        end
    end
      
    %u is in assigned, v is in unassigned
    if uindex <= length(setA)
        %u is in A
        setB = [setB; v];
    else
        %u is in B
        setA = [setA; v];
    end
    
    assigned = [setA; setB];
    unassigned = setdiff(unassigned, assigned);
    
     for i = 1:length(assigned)
        ku = assigned(i);
        G(ku, v) = 0;
        G(v, ku) = 0;
     end
    
    numassigned = length(assigned);
    
end

%if still some unassigned left, assign arbitrarily (to A, here)
% Basically, these have all zero-weight edges. i.e. zerocols
setA =[setA; unassigned'];
toc

nA = D(:, setA);
nB = D(:, setB);
nNa = length(setA);
nNb = length(setB);

nmu_a = coherence(nA);
nmu_b = coherence(nB);
nmu_larger  = max(nmu_a, nmu_b);
nmu_smaller = min(nmu_a, nmu_b);
    
nmu_m = mutualcoherence(nA,nB);
nmu_d = coherence(D);

nmu_d_alt = (0.5 * (nmu_larger*(k-2) + k*sqrt(nmu_larger^2 + nmu_m^2)))/(k-1);
nmu_d_eff = min(nmu_d, nmu_d_alt);

nk_bound1 = 2*(1 + nmu_larger)/(nmu_larger + 2*nmu_d_eff + sqrt(nmu_larger^2 + nmu_m^2));
nk_bound2 = (1+nmu_d_eff)/(2*nmu_d_eff);
nk_bound = max(nk_bound1, nk_bound2);

ndeltahat1 = 0.5 * (nmu_larger*(k-2) + k*sqrt(nmu_larger^2 + nmu_m^2));
ndeltahat2 = nmu_d * (k-1);
ndeltahat = min(ndeltahat1, ndeltahat2);
%ndeltahat =0;

valueset_optimized = [nNa, nNb, nmu_a, nmu_b, nmu_m, nmu_d, nmu_d_alt, nmu_d_eff, nk_bound, ndeltahat];

disp(paramlist);
disp(valueset_original);
disp(valueset_optimized);

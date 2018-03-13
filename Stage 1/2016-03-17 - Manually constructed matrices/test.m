
clear;clc;close all;
%perfectly sparse signals
%rng(10);

m = 80; %we're later taking m>= N/2
N = 100; 
nw = 20;

%m = 3;
%N = 5;
%nw = 3;


epsilon = 0;%1e-4; % for optimization threshold
epsilon_noise = 0;



%sensing matrix
%D = normrnd(1, 1, [m, N]);
%D = normc(D);
D1 = dctmtx(m);
D2 = eye(m);

D = [D1,D2];
D = D(:,randperm(size(D,2)));
D = D(:,1:N);


%% BOUNDS

%% 0. Direct
mu_d = coherence(D);
nx_bound0 = 0.5 * (1 + (1/mu_d));
delta0 = mu_d * (nw - 1);

%% 1. Arbitrary Direct Split dictionary
Na1 = round(N/2);
Nb1 = N - Na1;
D_A = D(:,1:Na1);
D_B = D(:,(Na1+1):N);

mu_a = coherence(D_A);
mu_b = coherence(D_B);
 
mu_m = mutualcoherence(D_A, D_B);
mu_larger = max(mu_a, mu_b);

nx_bound1 = 2*(1+mu_larger) / (mu_larger + (2*mu_d) + sqrt(mu_larger^2 + mu_m^2));
delta1 = 0.5 * (mu_larger*(nw-2) + nw*sqrt(mu_larger^2 + mu_m^2));

%% 2. Optimal place to split, no rearrangement

nx_bound2 = 0;
delta2 = Inf;
bestsplitpoint = 1; %1:1,2:N
for Na2 = 1:(N-1)
    Nb2 = N - Na2;
    D_A = D(:,1:Na2);
    D_B = D(:,(Na2+1):N);

    mu_a = coherence(D_A);
    mu_b = coherence(D_B);
    
    
    mu_m = mutualcoherence(D_A, D_B);
    mu_larger = max(mu_a, mu_b);
    
    candidatebound = 2*(1+mu_larger) / (mu_larger + (2*mu_d) + sqrt(mu_larger^2 + mu_m^2));
    candidatedelta = 0.5 * (mu_larger*(nw-2) + nw*sqrt(mu_larger^2 + mu_m^2));
    %disp(candidatebound);
    %if candidatebound > nx_bound2
    if candidatedelta < delta2    
        nx_bound2 = candidatebound;
        bestsplitpoint = Na2;
        delta2 = candidatedelta;
    end
end



%% 3. Optimally rearranged & split  (2-way)
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
setA =[setA; unassigned];

toc

D_A = D(:, setA);
D_B = D(:, setB);
Na3 = length(setA);
Nb3 = length(setB);

mu_a = coherence(D_A);
mu_b = coherence(D_B);
mu_m = mutualcoherence(D_A, D_B);
mu_larger = max(mu_a, mu_b);

nx_bound3 = 2*(1+mu_larger) / (mu_larger + (2*mu_d) + sqrt(mu_larger^2 + mu_m^2));
delta3 = 0.5 * (mu_larger*(nw-2) + nw*sqrt(mu_larger^2 + mu_m^2));




%% Display

disp('0: Dict as is')
disp(['   nx bound  = ' num2str(nx_bound0)]);
disp(['   delta_n_w = ' num2str(delta0)]);
disp(' ');

disp('1: dict split into halves');
disp(['   nx bound  = ' num2str(nx_bound1)]);
disp(['   delta_n_w = ' num2str(delta1)]);
disp(' ');

disp('2: dict split at optimal point, no shuffling');
disp(['   nx bound  = ' num2str(nx_bound2)]);
disp(['   delta_n_w = ' num2str(delta2)]);
disp(' ');

disp('3: optimal 2-way split');
disp(['   nx bound  = ' num2str(nx_bound3)]);
disp(['   delta_n_w = ' num2str(delta3)]);
disp(' ');



%% RECOVER x

% x_estimated = zeros(N,1);

l_reconstructionerrors = zeros(10);

ctr = 0;
for nw = 3:10
    ctr = ctr + 1;
    %signal
    x_actual = zeros(N,1); 
    x_actual(randsample(N, nw))= normrnd(0, 1, [nw, 1]);


    %measurement
    y = D * x_actual;

    noise = normrnd(0,1,[m,1]);
    noise = noise * epsilon_noise / norm(noise, 2);
    y_measured = y + noise;

    cvx_begin
        variables x_estimated(N)
        minimize(norm(x_estimated, 1))
        subject to
            norm(y_measured - (D*x_estimated)) <= epsilon 
            %norm(y_measured - (D*x_estimated) <= epsilon)
    cvx_end

    errors = (x_actual - x_estimated);
    reconstructionerror = sum(abs(errors));
    l_reconstructionerrors(ctr) = reconstructionerror;


    figure;
    plot(1:N, x_actual, '-r', 1:N, x_estimated, '-b');
    
    legend('actual', 'recovered');
    
end

l_reconstructionerrors(1:ctr)


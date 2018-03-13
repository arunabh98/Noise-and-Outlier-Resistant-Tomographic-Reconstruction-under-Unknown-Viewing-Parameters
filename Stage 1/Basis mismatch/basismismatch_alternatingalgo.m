%% PARAMETERS

clear; clc; close all;
N = 100; M = 60; sparse = 20;
r = 1e-1;
desirednormx = 100;
tolerance = 1e-4; % difference in consecutive estimates of x
epsilon = 1e-1; % for optimization threshold
epsilon_noise = 0;%epsilon;
maxiter = 100;



%% GENERATE DATA

%general
%x = normrnd(0,5,N,1);
%sparse
     x = zeros(N,1);
%compressible
%x = normrnd(0,0.1,N,1);
support = randsample(N,sparse);
x(support) = normrnd(0,5,sparse,1);
x = x * desirednormx / norm(x,2);

U = -(N-1)/2:(N-1)/2;
U = U';
deltaU = r*2*(rand(N,1) - 0.5);
Utilde = U + deltaU;
%deltaU = deltaU(1:M);
Delta = diag(deltaU);


%Assemble F on -(N-1)/2:(N-1)/2;
k=-(N-1)/2:(N-1)/2;
Psi=zeros(N,N);
Psitrue=zeros(N,N);
B=zeros(N,N);
for row = 1:N
    for col = 1:N
        Psitrue(row,col) = cos(2*pi * Utilde(row) * k(col)/N) * (1/sqrt(N));
        Psi(row,col) = cos(2*pi * U(row) * k(col)/N) * (1/sqrt(N));
        B(row,col) = -sin(2*pi * U(row) * k(col)/N) * (1/sqrt(N))*2*pi*k(col)/N;
    end
end

phi = normrnd(0,1,M,N);
yExact = phi*Psitrue*x;
yApprox = phi*(Psi + (B*Delta))*x

norm(yApprox-yExact)/norm(yExact)

noise = normrnd(0,1,[M,1]);
noise = noise * epsilon_noise / norm(noise, 2);
%y_measured = yExact + noise; %modeling error pre
y_measured = yApprox + noise; %no modeling error


%% RECOVERY

%alternating minimization

beta_estimated = zeros(N,1);

for iter = 1:50
    iter
    %find best x_estimated
    lambda=1;
    cvx_begin quiet
        variable x_estimated(N)
        minimize norm(x_estimated, 1) + lambda*norm(y_measured - (phi*(Psi + B*diag(beta_estimated))*x_estimated), 2)
    cvx_end

    %find best delta
    cvx_begin quiet
        variable beta_estimated(N)
        minimize norm(y_measured - (phi*(Psi+B*diag(beta_estimated))*x_estimated), 2)
    cvx_end
end


 
 delta_estimated = beta_estimated;
 delta_estimated(isnan(delta_estimated))=0;
 delta_estimated(delta_estimated>r) = r;
 delta_estimated(delta_estimated<-r) = -r;


%% RESULTS

%SIGNAL
figure();
plot(x, '-g');
hold on;
plot(x_estimated, '-r');
hold off;
legend('Actual', 'Estimated');

%DELTA U
figure();
plot(deltaU(support), '-b');
hold on;
plot(delta_estimated(support), '-r');
hold off;
legend('Actual', 'Estimated');


norm(x - x_estimated) / norm(x)
norm(deltaU(support) - delta_estimated(support)) / norm(deltaU(support))
norm(yApprox - (phi*(Psi+B*diag(delta_estimated))*x_estimated))/norm(y_measured)
norm(yApprox-yExact)/norm(yExact)


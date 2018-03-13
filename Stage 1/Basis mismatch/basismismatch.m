%% PARAMETERS

clear; clc; close all;
N = 100; M = 200; sparse = 10;
r = 1e-2;
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
yApprox = phi*(Psi + (B*Delta))*x;

norm(yApprox-yExact)/norm(yExact)

noise = normrnd(0,1,[M,1]);
noise = noise * epsilon_noise / norm(noise, 2);
%y_measured = yExact + noise; %modeling error pre
y_measured = yApprox + noise; %no modeling error


%% RECOVERY

%extended vector recov

W = [phi*Psi, phi*B];
%V = [x;Delta*x]
%y ~ W*V

lambda=10;

cvx_begin
    variable V(2*N)
    minimize norm(V, 1) + lambda*norm(y_measured - W*V, 2)
    %minimize norm(y_measured - W*V, 2)
cvx_end

%V=[x;Delta*x];

 
 x_estimated = V(1:N);
 delta_estimated = V(N+1:2*N) ./ x_estimated;
 delta_estimated(isnan(delta_estimated))=0;
 delta_estimated(delta_estimated>r) = r;
 delta_estimated(delta_estimated<-r) = -r;


%% RESULTS

signalerror = norm(x - x_estimated) / norm(x)
deltaerror = norm(deltaU(support) - delta_estimated(support)) / norm(deltaU(support))
norm(yApprox - W*V)/norm(y_measured)
taylorerror = norm(yApprox-yExact)/norm(yExact)


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

for i = 1:2*N
    W(:,i) = W(:,i)/norm(W(:,i));
end
Wcoh = max(max(abs(W'*W-eye(2*N))))

    
for i = 1:N
    phi(:,i) = phi(:,i)/norm(phi(:,i));
end
phicoh = max(max(abs(phi'*phi-eye(N))))



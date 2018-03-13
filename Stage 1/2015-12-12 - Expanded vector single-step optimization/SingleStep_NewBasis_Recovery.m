
clear; clc; close all;
rng(1);


%% PARAMETERS
m = 80;
n = 200;
kX = 10; % k-sparse signal
kBeta = 20;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
maxiter = 1000;


%% GENERATE DATA

A = normrnd(2, 2, [m, n]);
B = normrnd(2, 2, [m, n]);
U = dctmtx(n);

theta_actual = zeros(n,1);
theta_actual(randsample(n, kX)) = normrnd(0, 1, [kX, 1]);
beta_actual = zeros(n,1);
beta_actual(randsample(n, kBeta))= r*(2*rand(kBeta,1) - 1);

delta_actual = diag(beta_actual);
phi_actual = (A + (B * delta_actual));
y_actual = phi_actual * U * theta_actual;

%noise = zeros(m, 1);
noise = normrnd(0,0.5,[m,1]);
noise = noise * epsilon_noise / norm(noise, 2);

y_measured = y_actual + noise;


%% ESTIMATE


beta_estimated = zeros(n,1);
beta_estimated = r*(2*rand(n,1) - 1);
theta_estimated = zeros(n,1);
z_estimated = zeros(n,1);


rel_diff_theta = Inf;
R = r * ones(n,1);
iter = 0;

cvx_begin
    variables theta_estimated(n) z_estimated(n)
    minimize(norm([theta_estimated; z_estimated], 1))
    subject to
        norm(y_measured - [A*U B]*[theta_estimated; z_estimated], 2) <= epsilon
cvx_end

plot(1:n, theta_actual, '-r', 1:n, theta_estimated, '-b')









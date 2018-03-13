
clear; clc; 
rng(5);

%% PARAMETERS
n = 200;
r = 0.1;
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-2; % for optimization threshold
epsilon_noise = epsilon;
numruns = 10;
mlist = 10:10:150;
klist = 10:5:60;

%% INITIALIZE

x_error_averaged = zeros(numel(mlist), numel(klist));

for run = 1:numruns;

    disp('====================================');
    disp(['              run ' num2str(run) '                 ']);
    x_error = zeros(numel(mlist), numel(klist));

    for mcounter = 1:numel(mlist)

        m = mlist(mcounter) % number of measurements

        %% GENERATE DATA

        A = normrnd(2, 5, [m, n]);
        B = normrnd(2, 5, [m, n]);

        beta_actual = r*(2*rand(n,1) - 1); 
        delta_actual = diag(beta_actual);
        phi_actual = (A + (B * delta_actual));

        x_actual = zeros(n,1);

        for kcounter = 1:numel(klist)% k-sparse signal
            k = klist(kcounter)

            x_actual_k = x_actual;
            x_actual_k(randsample(n, k)) = normrnd(0, 1, [k, 1]);

            y_actual = phi_actual * x_actual_k;

            noise = zeros(m, 1);
            noise = normrnd(0, 1, [m,1]);
            noise = noise * epsilon_noise / norm(noise, 2);

            y_measured = y_actual + noise;


            %% ESTIMATE

            beta_estimated = zeros(n,1);
            x_plus = zeros(n,1);
            x_minus = zeros(n,1);
            p = zeros(n,1);
            M = [A -A B];
            R = r * ones(n,1);

            cvx_begin
                variables x_plus(n) x_minus(n) p(n)
                minimize(sum(x_plus + x_minus))

                subject to
                    norm(y_measured - M*[x_plus; x_minus; p], 2) <= epsilon

                    x_plus >= 0
                    x_minus >= 0
                    -r*(x_plus + x_minus) <= p <= r*(x_plus + x_minus) 
            cvx_end

            assert(all(x_plus .* x_minus < 1e-8), 'CHECK FAIL: x+ and x- overlap')
            x_estimated = x_plus - x_minus;
            beta_estimated(x_estimated > 1e-8) = p(x_estimated > 1e-8) ./ x_estimated(x_estimated > 1e-8); % Avoid div by 0 locations

            x_error(mcounter, kcounter) = norm(x_estimated - x_actual_k, 1);



        end % k loop

    end % m loop

    x_error_averaged = x_error_averaged  + x_error;
    
end % run loop

x_error_averaged = x_error_averaged / numruns;    

%% EVALUATE RESULTS

% figure;
% plot(1:n, x_actual_k, '-r', 1:n, x_estimated, '-b');
% legend('actual', 'recovered');
% 
% norm(x_actual_k - x_estimated, 1)    % AA-
% 
%     

beep; pause(1); beep; pause(1); beep; 
save('x_error_averaged.mat', 'x_error_averaged');

image(klist, mlist, repmat(x_error_averaged / max(max(x_error_averaged)),[1,1,3]))
xlabel('k (number of non-zero entries)');
ylabel('m (number of measurements)');
title('L1 norm of error vector')



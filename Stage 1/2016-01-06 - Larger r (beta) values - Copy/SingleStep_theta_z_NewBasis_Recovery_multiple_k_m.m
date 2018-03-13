
clear; clc; 
rng(5);

%% PARAMETERS

n = 200;
                 
tolerance = 1e-6; % difference in consecutive estimates of x
epsilon = 1e-1; % for optimization threshold
epsilon_noise = epsilon;

numruns = 5;
% mlist = 10:20:150;
% kThetalist = 10:10:60;
% kBetalist = 15:15:90;
% rlist = 0.5:0.5:4;

mlist = 100;
kThetalist = 20;
kBetalist = 20;
rlist = 0.5


%% INITIALIZE


for r = rlist
    R = r * ones(n,1);
    theta_error_averaged = zeros(numel(mlist), numel(kThetalist), numel(kBetalist));

    tic
    for run = 1:numruns;

        disp('====================================');
        disp(['              run ' num2str(run) '                 ']);
        theta_error = zeros(numel(mlist), numel(kThetalist), numel(kBetalist));

        for mcounter = 1:numel(mlist)

            m = mlist(mcounter) % number of measurements

            %% GENERATE DATA

            A = normrnd(2, 2, [m, n]);
            B = normrnd(2, 2, [m, n]);
            U = dctmtx(n);

            theta_actual = zeros(n,1);
            beta_actual = zeros(n,1);

            noise = normrnd(0,1,[m,1]);
            noise = noise * epsilon_noise / norm(noise, 2);


            for kThetacounter = 1:numel(kThetalist)% signal sparsity
                kTheta = kThetalist(kThetacounter)

                theta_actual_k = theta_actual;
                theta_actual_k(randsample(n, kTheta)) = normrnd(0, 1, [kTheta, 1]);

                for kBetacounter = 1:numel(kBetalist)% k-sparse signal  
                    [r run mcounter kThetacounter kBetacounter]

                    kBeta = kBetalist(kBetacounter)
                    beta_actual_k = beta_actual;

                    beta_actual_k(randsample(n, kBeta))= r*(2*rand(kBeta,1) - 1);

                    delta_actual_k = diag(beta_actual_k);
                    phi_actual_k = (A + (B * delta_actual_k));

                    y_actual = phi_actual_k * U * theta_actual_k;
                    y_measured = y_actual + noise;


                    %% ESTIMATE

                    beta_estimated = zeros(n,1);
                    beta_estimated = r*(2*rand(n,1) - 1);
                    theta_estimated = zeros(n,1);
                    z_estimated = zeros(n,1);

                    cvx_begin
                        variables theta_estimated(n) z_estimated(n)
                        minimize(norm([theta_estimated; z_estimated], 1))
                        subject to
                            norm(y_measured - [A*U B]*[theta_estimated; z_estimated], 2) <= epsilon
                    cvx_end

                    theta_error(mcounter, kThetacounter, kBetacounter) = norm(theta_estimated - theta_actual_k, 1);
                end % kBeta loop

            end % kTheta loop

        end % m loop

        theta_error_averaged = theta_error_averaged  + theta_error;

    end % run loop
    toc
    theta_error_averaged = theta_error_averaged / numruns;    


    save(['theta_error_averaged_' num2str(r) '.mat'], 'theta_error_averaged');

    %% PLOT & EVALUATE RESULTS

    f1 = figure;
    for i = 1:numel(mlist) 
        subplot(2,4,i)
        btmat = squeeze(theta_error_averaged(i,:,:))
        image(kBetalist, kThetalist, repmat(btmat / max(max(btmat)),[1,1,3]))
        ylabel('kTheta');
        xlabel('kBeta');
        title(['m = ' num2str(mlist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized separately for each individual plot');
    f2 = figure;
    for i = 1:numel(mlist) 
        subplot(2,4,i)
        btmat = squeeze(theta_error_averaged(i,:,:))
        image(kBetalist, kThetalist, repmat(btmat / max(max(max(theta_error_averaged))),[1,1,3]));
        ylabel('kTheta');
        xlabel('kBeta');
        title(['m = ' num2str(mlist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized over all m');    

    %%%%%%%%%%%%%%%

    f3 = figure;
    for i = 1:numel(kThetalist) 
        subplot(2,3,i)
        mbmat = squeeze(theta_error_averaged(:,i,:))
        image(kBetalist, mlist, repmat(mbmat / max(max(mbmat)),[1,1,3]))
        ylabel('m');
        xlabel('kBeta');
        title(['kTheta = ' num2str(kThetalist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized separately for each individual plot');

    f4 = figure;
    for i = 1:numel(kThetalist) 
        subplot(2,3,i)
        mbmat = squeeze(theta_error_averaged(:,i,:))
        max(max(mbmat))
        image(kBetalist, mlist, repmat(mbmat /  max(max(max(theta_error_averaged))),[1,1,3]))
        ylabel('m');
        xlabel('kBeta');
        title(['kTheta = ' num2str(kThetalist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized over all kThetas');


    %%%%%%%%%%%%%%%%%%

    f5 = figure;
    for i = 1:numel(kBetalist) 
        subplot(2,3,i)
        mtmat = squeeze(theta_error_averaged(:,:,i))
        image(kThetalist, mlist, repmat(mtmat / max(max(mtmat)),[1,1,3]))
        ylabel('m');
        xlabel('kTheta');
        title(['kBeta = ' num2str(kBetalist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized separately for each individual plot');

    f6 = figure;
    for i = 1:numel(kBetalist) 
        subplot(2,3,i)
        mtmat = squeeze(theta_error_averaged(:,:,i))
        image(kThetalist, mlist, repmat(mtmat / max(max(max(theta_error_averaged))),[1,1,3]))
        ylabel('m');
        xlabel('kTheta');
        title(['kBeta = ' num2str(kBetalist(i))]);
    end
    suptitle('L1 norm of error vector. Normalized over all kBetas');

    
    %% SAVE EVERYTHING
    save(['env_' num2str(r) '.mat']);
    print(f1, ['r=' num2str(r) '_f1_varying_m.png'],'-dpng');
    print(f2, ['r=' num2str(r) '_f2_varying_m2.png'],'-dpng');
    print(f3, ['r=' num2str(r) '_f3_varying_kTheta.png'],'-dpng');
    print(f4, ['r=' num2str(r) '_f4_varying_kTheta2.png'],'-dpng');
    print(f5, ['r=' num2str(r) '_f5_varying_kBeta.png'],'-dpng');
    print(f6, ['r=' num2str(r) '_f6_varying_kBeta2.png'],'-dpng');
    
    close all;
    
end 

figure;
plot(1:n, theta_actual_k, '-r', 1:n, theta_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('theta actual', 'theta recovered');

beta_estimated = z_estimated ./ (U * theta_estimated);
figure;
plot(1:n, beta_actual_k, '-r', 1:n, beta_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('beta actual', 'beta recovered');

figure;
plot(1:n, beta_actual_k .* (U * theta_actual_k), '-r', 1:n, z_estimated, '-b');
title(['m = ' num2str(m) ' | kTheta = ' num2str(kTheta) ' | kBeta = ' num2str(kBeta)]);
legend('z actual', 'z recovered');



%beep; pause(1); beep; pause(1); beep; 



%function numcorrect = ARP(imagefile, sigmaNoiseFraction, numkeep, numstarts, Ord, randseed)


    % Reconstruction from projection data from unknown angles
    % In this script we'll implement angle recovery. From there, reconstruction
    % is quite directly possible

    close all; clear; clc;
    randseed = 1;
    rng(randseed);

    tic
    %% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
    thetas = (0:179); % Projections to sample from
    numkeep = 30; % Number of angles to recompute from
    sigmaNoiseFraction = 0.00;

    %% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

    %%%% Define image I %%%%
    %I = phantom(200);
    imagefile = '../images/200px-mickey.jpg';
    I = rgb2gray(imread(imagefile));
    %  figure();
    %imshow(I,'InitialMagnification','fit');


    %%%% Calculate projections %%%%
    [P, svector] = radon(I, thetas); % P(s, theta)

    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;

    % Sample Uniformly
    keep = sample(P, numkeep); %returns indices
    numkeep = length(keep); %to adjust for rounding off errors not balancing out

    thetaskeep = thetas(keep);
    Pgiven = P(:, keep);
    Pnonoise = Pgiven;

    %expectedanswer = (thetaskeep - thetaskeep(1))';
    expectedanswer = thetaskeep';
    expectedanswer  = expectedanswer * mode(sign(expectedanswer));

    %%%% Add noise
    ref = std(Pgiven(:));
    sigmaNoise = sigmaNoiseFraction * ref;
    noise = normrnd(0, sigmaNoise, size(Pgiven));
    Pgiven = Pgiven + noise;

   
    
    %% %%%%%%%%%%% B. DENOISING %%%%%%%%
    if sigmaNoiseFraction > 0
        
        disp('Before denoising');
        Pnoisy = Pgiven;
        norm(Pnonoise - Pgiven, 'fro')
        psnr_noisy = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - Pnoisy ) .^ 2))) )

        P_denoised = denoise(Pgiven, sigmaNoise, 50, 200);
        Pgiven = P_denoised;

        disp('After denoising')
        norm(Pnonoise - P_denoised, 'fro')
        %norm(Pnonoise - P_denoised2, 'fro')
        figure()
        plot(Pnonoise(:,1));
        hold on;
        plot(P_denoised(:,1), '-r');
        plot(Pnoisy(:,1), '-g');
        hold off;

        psnr_denoised = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - P_denoised ) .^ 2))) )

    end
    %% %%%%%% B2. Final fixes on P, M %%%%%%%%%%%%%%%
    Pgiven = max(0, Pgiven);

    %% %%%%%%%%%%%%%%% C. RECOVER ANGLES %%%%%%%%%%%%%%%%%

    % Initialize angles randomly
    thetasestimated = randi([1 179], numkeep, 1);
    %thetasestimated(1) = 0;
    thetasfirstestimate = thetasestimated;

    
    
    %reconstruct with actual
    I_reconstructed_original = iradon(Pgiven, thetaskeep);
    %figure();
    %imshow(I_reconstructed_original/max(max(I_reconstructed_original)), 'InitialMagnification','fit');
    
    %reconstruct with current estimate 
    I_reconstructed = iradon(Pgiven, thetasestimated);
    %figure();
    %imshow(I_reconstructed / max(max(I_reconstructed)), 'InitialMagnification','fit');
    
    
    
    D_original_angles = dct2(I_reconstructed_original);
    norm(D_original_angles, 1);
    sparsity(D_original_angles)
        
    D_current = dct2(I_reconstructed);
    norm(D_current, 1);
    sparsity(D_current)
    
    D_actual_image = dct2(I);
    norm(D_actual_image, 1);
    sparsity(D_actual_image)
    
    D_zeros = dct2(iradon(Pgiven, zeros(size(thetaskeep))));
    norm(D_zeros, 1);
    sparsity(D_zeros)
    
    
    %F1 = fft2(I_reconstructed_original);
    %F2 = fft2(I_reconstructed);
       
   
    I_current = iradon(Pgiven, thetasestimated);
    alpha_current = dct2(I_current);
    
    E = norm(alpha_current, 1);
    %E = sparsity(alpha_current);
    E_previter = Inf;
    
    numiter = 30;
    E_values = zeros(numiter, 1);
    
    ctr = 0;
    for iteration = 1:numiter
        ctr = ctr + 1;
        disp([num2str(iteration)]);
        for i = 1:numkeep
            besttheta = thetasestimated(i);
            bestE = E;
            for t = [-179:180]
                thetas_iter = thetasestimated;
                thetas_iter(i) = t;

                I_iter = iradon(Pgiven, thetas_iter);
                alpha_iter = dct2(I_iter);
                E_t = norm(alpha_iter, 1);

                if E_t < bestE
                    besttheta = t;
                    bestE = E_t;
                end
            end

            thetasestimated(i) = besttheta;
            I_current = iradon(Pgiven, thetasestimated);
            alpha_current = dct2(I_current);
            E = norm(alpha_current, 1);
            disp(['theta', num2str(i), ' | norm1 = ', num2str(E)]);

        end
        
        E_iter = E;
        E_values(ctr) = E_iter;
        delta = E_previter - E_iter;
        if delta < 0.0001
            break
        end
        disp(['--- End iter', num2str(ctr), '---'])
        E_previter = E_iter;
    end
    
    
    
    
    

    I_reconstructed = iradon(Pgiven, thetasestimated);
    alpha_reconstructed  = dct2(I_reconstructed );
    norm(alpha_reconstructed, 1)
    sparsity(alpha_reconstructed)
    
    figure();
    imshow(I_reconstructed / max(max(I_reconstructed)), 'InitialMagnification','fit');
    
    
    
    
    
    
    %end
    
    

    
    
    %% %%%%%%%%%%%%%%%%% IMAGE RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%

% %     %Using ORIGINAL angles
% %     I_reconstructed_original = iradon(Pgiven, thetaskeep);
% %     figure();
% %     imshow(I_reconstructed_original/max(max(I_reconstructed_original)), 'InitialMagnification','fit');
% % 
% % 
% %     %Using ALL RECONSTRUCTED angles (including incorrect ones)
% %     I_reconstructed = iradon(Pgiven, thetasestimated);
% %     figure();
% %     imshow(I_reconstructed/max(max(I_reconstructed)),'InitialMagnification','fit');
% % 
% %     % a(logical([1,0,1,1,0,0,0,1,0,0])) = []
% %     beep; pause(1); beep; pause(1); beep;

%end

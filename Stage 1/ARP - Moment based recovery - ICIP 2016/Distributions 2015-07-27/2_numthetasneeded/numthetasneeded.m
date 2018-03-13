%% How many angles do we need to get a decent reproduction
clear;clc;

resolution = 1;
numkeeprange = 10:10:(179/resolution);
thetas = (0:resolution:179); % Projections to sample from

I = phantom(200); imagename= 'phantom';
% I = ibeam(); imagename = 'ibeam';
% I = rgb2gray(imread('../images/256px-Earthrise.jpg')); imagename = 'earthrise';
% I = rgb2gray(imread('../images/200px-beach.jpg')); imagename = 'beach';
% I= rgb2gray(imread('../images/200px-woody.jpg')); imagename = 'woody';


[P, svector] = radon(I, thetas); % P(s, theta)
smax = max(abs(svector));
svector = svector / smax;

% Adding noise
sigmaNoiseFraction = 0.05;

ref = std(P(:));
sigmaNoise = sigmaNoiseFraction * ref;
noise = normrnd(0, sigmaNoise, size(P));
P = P + noise;
P = max(0, P);

%outdir = ['output\noise5pc\' 'resolution' num2str(resolution) '\'];
outdir = 'output\alter1in25thetas_5pcnoise\';

imwrite(I, [outdir imagename '_0original.bmp']);
for numkeep = numkeeprange
    if mod(numkeep, (1/resolution)) == 0
        numkeep
    end
    
    keep = randperm(size(P,2), numkeep);
    keep = sort(keep);
    thetaskeep = thetas(keep);
    Pgiven = P(:, keep);
    
    % Altering thetas
    numalter = round(numkeep/25);
    thetasalter = randperm(numkeep, numalter);
    thetaskeep(thetasalter) = randi(180, numalter, 1); % assign an arbitrary value

    % Reconstruct
    I_reconstructed = iradon(Pgiven, thetaskeep);
    I_reconstructed = I_reconstructed/max(max(I_reconstructed));
    
    %imshow(I_reconstructed,'InitialMagnification','fit');
    imwrite(I_reconstructed, [outdir imagename '_' num2str(numkeep) '.bmp']);
end


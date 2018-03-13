function P_denoised = denoise2(Pgiven, sigmaNoise, patchsize, L)

%PCA based denoising 


%Pgiven, sigmaNoise

% % Pgiven = repmat(10, 50, 10);
% % Pnonoise = Pgiven;
% % noise = normrnd(0, 0.1, size(Pgiven));
% % sigmaNoise = std(noise(:));
% % Pgiven = Pgiven + noise;


% patchsize = 10;
% L = 5; % number of similar patches to consider

patchlist = zeros(patchsize, size(Pgiven, 2) * (size(Pgiven, 1) - patchsize + 1) );
estimate = zeros(size(Pgiven));
countestimates = zeros(size(Pgiven));



ctr = 0;
for theta = 1:size(Pgiven, 2)
    for patchstart = 1:(size(Pgiven, 1) - patchsize + 1)
        ctr = ctr + 1;
        patch = Pgiven(patchstart:(patchstart+patchsize-1), theta);
        patchlist(:, ctr) = patch;
    end
end


patchlist_denoised = zeros(size(patchlist));

for p = 1:ctr
    if mod(p,100) == 0
        fprintf(' %d', p);
    end
    qref = patchlist(:,p); %column vector
    differences = patchlist - repmat(qref, 1, ctr);
    distances = sqrt(sum(differences .^ 2));
    
    [~, distanceorder] = sort(distances);
%     nearestneightbours = distanceorder(1:(L+1)); % nearest L+1 neighbours
%     nearestneightbours = setdiff(nearestneightbours , p); %excluding p itself
%     nearestneightbours = nearestneightbours(1:L); %in the weird case where p wasn't picked at all

    nearestneightbours = distanceorder(1:L); % nearest L neighbours
    xref = patchlist(:, nearestneightbours); % qref is a part of this
    
    [evec, eval] = eig(xref * xref');
    
    
    coeffref = evec' * qref;
    coeffs = evec' * xref;
    % so that now:     evec * coeffs(i) ~= xref(:, i)
    

    abarsq = (sum(coeffs .^ 2, 2) / L) - (sigmaNoise^2);
    abarsq = max(0, abarsq);
    b_ref = coeffref ./ (1 + (sigmaNoise ^ 2) ./ abarsq);
    b_all = coeffs   ./ (1 + (sigmaNoise ^ 2) ./ repmat(abarsq, 1, L));
    
    %Denoise only reference
    qdenoised = evec * b_ref;
    
    %Denoise everything in eigenspace
    P = evec * b_all;
    patchespertheta = size(Pgiven, 1) - patchsize + 1;
    
    for n = 1:L
        p2 = nearestneightbours(n);
        theta  = ceil(p2 / patchespertheta);
        start_s = 1 + mod(p2 - 1, patchespertheta);
        for i = 0:patchsize - 1
            s = start_s + i;
            estimate(s, theta) = ((estimate(s, theta) * countestimates(s, theta)) + P(i + 1, n)) / (countestimates(s, theta) + 1);
            countestimates(s, theta) = countestimates(s, theta) + 1;
        end
    end
    
    
    
    patchlist_denoised(:, p) = qdenoised;
    
        
end
disp(' ')
disp([num2str(ctr) ' total patches denoised'])


%Reassemble P from patches


% ctr = 0;
% P_denoised = zeros(size(Pgiven));
% numestimates = zeros(size(P_denoised));
% 
% for theta = 1:size(P_denoised, 2)
%     for patchstart = 1:(size(P_denoised, 1) - patchsize + 1)
%         ctr = ctr + 1;
%         patch = patchlist_denoised(:, ctr);
%         P_denoised(patchstart:(patchstart+patchsize-1), theta) =  P_denoised(patchstart:(patchstart+patchsize-1), theta) + patch;
%         numestimates(patchstart:(patchstart+patchsize-1), theta) = numestimates(patchstart:(patchstart+patchsize-1), theta) + 1;
%     end
% end
% 
% P_denoised = P_denoised ./ numestimates;
P_denoised = estimate;





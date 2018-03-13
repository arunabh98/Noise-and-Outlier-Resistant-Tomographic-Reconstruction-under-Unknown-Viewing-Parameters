%denoising

%Pgiven, sigma
% Pgiven = repmat(10, 287, 31);
% Pnonoise=0
% Pnonoise = Pgiven;
% noise = normrnd(0, 1, size(Pgiven));
% Pgiven = Pgiven + noise;


patchsize = 25;
L = 50; % number of similar patches to consider

patchlist = zeros(patchsize, size(Pgiven, 2) * (size(Pgiven, 1) - patchsize + 1) );

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
    
    [distancessorted, distanceorder] = sort(distances);
%     nearestneightbours = distanceorder(1:(L+1)); % nearest L+1 neighbours
%     nearestneightbours = setdiff(nearestneightbours , p); %excluding p itself
%     nearestneightbours = nearestneightbours(1:L); %in the weitrd case where p wasn't picked at all

    nearestneightbours = distanceorder(1:L); % nearest L neighbours
    xref = patchlist(:, nearestneightbours); % SHOULD qref BE A PART OF THIS?
    
    [evec, eval] = eig(xref * xref');
    
    coeffs = zeros(patchsize, L);
    coeffref = evec' * qref;
    
    for i = 1:L
        coeffs(:, i) = evec' * xref(:, i);
        % so that now:     evec * coeffs(i) ~= xref(:, i)
    end

    abarsq = (sum(coeffs .^ 2, 2) / L) - (sigmaNoise^2);
    abarsq = max(0, abarsq);
    bref = coeffref ./ (1 + (sigmaNoise ^ 2)/abarsq)';
    
    qdenoised = evec * bref;
    patchlist_denoised(:, p) = qdenoised;
        
end
disp(' ')


%Reassemble P from patches


ctr = 0;
P_denoised = zeros(size(Pgiven));
numestimates = zeros(size(P_denoised));

for theta = 1:size(P_denoised, 2)
    for patchstart = 1:(size(P_denoised, 1) - patchsize + 1)
        ctr = ctr + 1;
        patch = patchlist(:, ctr);
        P_denoised(patchstart:(patchstart+patchsize-1), theta) =  P_denoised(patchstart:(patchstart+patchsize-1), theta) + patch;
        numestimates(patchstart:(patchstart+patchsize-1), theta) = numestimates(patchstart:(patchstart+patchsize-1), theta) + 1;
    end
end

P_denoised = P_denoised ./ numestimates;



norm(Pnonoise - Pgiven)
norm(Pnonoise - P_denoised)



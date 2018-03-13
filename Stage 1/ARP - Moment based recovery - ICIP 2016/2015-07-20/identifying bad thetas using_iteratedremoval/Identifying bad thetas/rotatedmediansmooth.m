function smoothed= rotatedmediansmooth(curve, thetareference_index, n)
    smoothed = smooth(curve, 'sgolay', 3);
    if nargin == 2
        n = 3;
    end
    numkeep = length(curve);
    
    curverotated = [curve(thetareference_index:numkeep) ; curve(1:(thetareference_index-1)) ];
    medianfilterrotated = medfilt1(curverotated, n);
    %reverse-rotate
    smoothed = [medianfilterrotated((numkeep - thetareference_index + 2):numkeep); ...
                    medianfilterrotated(1:(numkeep - thetareference_index + 1))];

    
                
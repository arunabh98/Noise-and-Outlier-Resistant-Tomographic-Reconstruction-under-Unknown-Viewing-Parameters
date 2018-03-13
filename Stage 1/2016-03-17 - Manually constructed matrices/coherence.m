function mu = coherence(M)
    M = normc(M);
    dotproducts = M'*M;
    dotproducts(logical(eye(size(dotproducts)))) = 0;
    mu = max(max(abs(dotproducts)));
end

    

function mu_m = mutualcoherence(A,B)
    A = normc(A);
    B = normc(B);
    mutualdotproducts = A'*B;
    mu_m = max(max(abs(mutualdotproducts)));
end

    

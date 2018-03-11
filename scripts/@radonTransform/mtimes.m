function res = mtimes(A, x)
D = dctmtx(A.height);
if A.adjoint == 0 %A*x
    x = reshape(x, [A.height, A.width]);
    x = D'*x;
    res = radon(x, A.angles);
    res = res(:);
else %At*x
    x = reshape(x, [A.projection_length, size(A.angles, 2)]);
    res = 2*iradon(x, A.angles, 'linear', 'none', A.output_size);
    res = D*res;
    res = res(:);
end
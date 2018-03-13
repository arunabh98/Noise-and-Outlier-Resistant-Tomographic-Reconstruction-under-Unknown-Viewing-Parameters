L_list = [10,100];
patchsize_list = [10:40:200];

L_mat = repmat(L_list, length(patchsize_list), 1);
patchsize_mat = repmat(patchsize_list, length(L_list), 1)';

z = zeros(size(patchsize_mat));

for i = 1:length(patchsize_list)
    for j = 1:length(L_list)
        L = L_mat(i,j)
        patchsize = patchsize_mat (i,j)
        P_denoised = denoise2(Pgiven, sigmaNoise, patchsize, L);
        z(i,j) = norm(Pnonoise - P_denoised, 'fro');
        
    end
end

figure();
surf(L_mat, patchsize_mat, z);
xlabel('L - #nearestneighbours');
ylabel('patchsize');
zlabel('residualnoise');

print('residualnoise2', '-djpeg')


plot(Pnonoise(:,1))
line(P_denoised(:,1))
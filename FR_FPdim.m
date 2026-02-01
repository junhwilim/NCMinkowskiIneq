function [d]=FR_FPdim(N)
% returns the Frobenius-Perron dimension 
% of the basis elements of the fusion ring N.
[r,s] = size(N); %r=1
d = zeros(1,s);
for i = 1:s
    d(1,i) = norm(N{i}); 
    %By Frobenius Reciprocity the dual and matrix adjoint coincide. 
    %Since the fusion matrices have the common PF eigenvector,
    %  PF eigenvalue = spectral radius = norm. (by C*-identity) 
end
end


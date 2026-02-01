function [M,d]=FR_fusion_matrix(N)
% N{i}(j,k)=N_{ij}^k (The input is the transpose of the fusion matrices.)
[r,s] = size(N); %r=1
M = cell(1,s);
for i=1:s
    M{i}=(N{i}).'; % N{i}(k,j)=N_{ij}^k 
end
[d]=FR_FPdim(N);
end
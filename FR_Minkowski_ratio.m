function [MinRat] = FR_Minkowski_ratio( N, d, s, i, j, k, ap )
% Computes RHS/LHS of the Minkowski inequality
% N: fusion ring, d: Frobenius-Perron dimension, s: rank
% 1<= i,j,k <= s: index of the basis element,
% ap: (input vector a) + (reciprocal of the input exponent p)
    a = ap(1,1:s); %input vector
    p = ap(s+1)^(-1); %input exponent
%% RHS
    RHS=0;
    for m = 1:s
        inside=0;
        for n = 1:s
            summand = ( N{i}(m,n) * (a(1,n)^p) * d(1,n) / d(1,m) ); %N[i][m,n]=N_{i,m}^n
            inside = inside + summand;
        end
        RHS = RHS + ( (inside^(1/p)) * N{j}(k,m) * d(1,m) / d(1,j) ); %N[j][k,m]=N_{j,k}^m
    end
    if RHS < 10^(-12)
        MinRat=0;
        return
    end
%% LHS
    LHS=0;
    for l = 1:s
        inside=0;
        for n = 1:s
            inside = inside + ( N{l}(k,n) * a(1,n) * d(1,n) / d(1,l) );  %N[l][k,n]=N_{l,k}^n
        end
        LHS = LHS + ( inside^p * N{i}(j,l) * d(1,l) / d(1,j) );  %N[i][j,l]=N_{i,j}^l
    end
    LHS = LHS^(1/p);
    if LHS < 10^(-12)
        fprintf("LHS too small\n")
        MinRat=0;    
        return
    end
%% RHS/LHS
    MinRat = RHS/LHS;
end

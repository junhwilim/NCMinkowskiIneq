function [Delx]=FR_Del_x( N, d, s, i, j, k, ap, Delta_x )
% Computes the gradient of RHS/LHS of the Minkowski inequality.
    Delx = zeros(1,s+1);
    for n = 1:s+1
        b1p = ap;
        b1p(1,n) = b1p(1,n) - Delta_x;
        if b1p(1,n) < 0
            fprintf("b1p(1,n)=%f", b1p(1,n));
        end
        MinRat1 = FR_Minkowski_ratio( N, d, s, i, j, k, b1p );
        
        b2p = ap;
        b2p(1,n) = b2p(1,n) + Delta_x;
        MinRat2 = FR_Minkowski_ratio( N, d, s, i, j, k, b2p );
        
        if MinRat1 == 0 || MinRat2 == 0
            Delx=0;
            return
        end
        partial_n = ( MinRat2 - MinRat1 ) / ( 2 * Delta_x );
        Delx(1,n) = partial_n;
    end

end

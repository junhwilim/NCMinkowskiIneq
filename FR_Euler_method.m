function [is_unitary,V] = FR_Euler_method( N, Delta_x, Delta_t, steps )
% Implements gradient descent to RHS/LHS of Minkowski inequality
% N: fusion ring, Delta_x: increment used for approximation of gradient
% Delta_t: step size, steps: # of steps
    V=cell(0);
    cell_index=0;
    is_unitary = 1; % will turn 0 if RHS/LHS becomes <1
    [r,s] = size(N); %r=1, s=rank
    pd = makedist('Uniform','lower',0,'upper',1); %uniform distribution on [0,1]
   
    [d]=FR_FPdim(N);
       
    for i = 1:s
        for j = 1:s
            for k = 1:s
                validity = false;
                while ~validity
                    ap = random(pd,[1,s+1]); %randomly generated vector in [0,1]^(s+1)
                    %first s entries: input vector a
                    %last entry: reciprocal of the input exponent p
                    validity = input_validity( ap, Delta_x, s ); %imposes restriction on the input
                end
                    
                for l = 1:steps
                    Delx = FR_Del_x( N, d, s, i, j, k, ap, Delta_x ); %gradient of RHS/LHS
                    if Delx == 0
                        %fprintf("Delx==0\n")
                        break
                    end
                    ap1 = ap - Delta_t * Delx ; %updated ap
                    if ~(input_validity( ap1, Delta_x, s ))
                        %fprintf("input_validity( ap1, Delta_x, s )==false\n")
                        break
                    end
                    %fprintf("i=%d j=%d k=%d\n",i,j,k);
                    ap = ap1;
                end
                MinRat = round( FR_Minkowski_ratio( N, d, s, i, j, k, ap ), 14 ); %RHS/LHS
                if MinRat < 1 && MinRat > 0 && ap(1,s+1) > 0.1 %i.e. fusion ring N excluded
                    is_unitary=0;
                    fprintf("i=%d, j=%d, k=%d, a=[", i, j, k);
                    for a_index =1:s-1
                        fprintf("%f,",ap(1,a_index))
                    end
                    fprintf("%f], p=%f, RHS/LHS=%f\n", ap(1,s), ap(1,s+1)^(-1), MinRat );
                    cell_index = cell_index+1;
                    V{cell_index,1}={i,j,k,ap,MinRat};
                end
            end
        end
    end
end

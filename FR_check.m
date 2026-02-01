function [validity]=FR_check(N)
% Checks whether N is a valid fusion ring.
% N{i}(j,k)=N_{ij}^k (The input is the transpose of the fusion matrices.)
validity=1;
[r,s] = size(N); %r=1
% converts the fusion ring into the set of fusion matrices
% i.e. M{i}(k,j)=N_{ij}^k
[M]=FR_fusion_matrix(N);

%% Size and Positivity
% Checking the size of the matrices
for i=1:s
    [m,n]=size(M{i});
    if m ~= s || n~=s
        fprintf("N{%d} is not a %d-by-%d matrix.\n",i,s,s);
        validity=0;
    end
end

% Checking the positivity of the entries
for i=1:s
    for j=1:s
        for k=1:s
            if M{i}(k,j)<0
                fprintf("N{%d}(%d,%d) is negative.\n",i,j,k);
                validity=0;
            end
        end
    end
end

%% Identity
% Checking N{1} is the identity
% -> N_{1j}^i = \delta_{ij} & 1^*=1
error = norm( M{1} - eye(s) );
if error>0
    fprintf("The first matrix is not the identity.\n");
    validity=0;
end

% Checking N_{i1}^j= \delta_{ij}
for i=1:s
    error = M{i}(:,1);
    error(i,1) = error(i,1)-1;
    if norm(error)>0
        fprintf("N{%d}(1,:) is not e_%d.\n",i,i);
        validity=0;
    end
end

%% Dual
% Checking N_{i*}^1=e_{i^*} for a unique i^*
dual = zeros(1,s);
for i=1:s
    for j=1:s
        if M{i}(1,j) > 0
            break
        end
    end
    error=M{i}(1,:);
    error(1,j) = error(1,j) - 1;
    if norm(error) > 0
        fprintf("N{%d}(:,1) is not e_{i^*} for any i^*.\n",i);
        validity=0;
    else
        dual(1,i) = j;
    end
end

% Checking the involutivity of the dual
% -> N_{i^*i}^1=1 && N_{i*}^1 are linearly independent
for i=1:s
    if dual(1,dual(1,i)) ~= i
        fprintf("dual(%d)=%d, dual(dual(%d)))=%d.\n",i,dual(1,i),i,dual(dual(1,i)));
        validity=0;
    end
end

%% Associativity
% Checking (N_{i*}^*)(N_{j*}^*)= \sum_k N_{ij}^k(N_{k*}^*)
% -> Associativity 
% (by the associativity of matrix multiplication & linear independence)
for i=1:s
    for j=1:s
        difference = M{i}*M{j};
        for k=1:s
            difference = difference - M{i}(k,j)*M{k};
        end
        error=norm(difference);
        if error>0
            fprintf("Error in fusion coefficients: i=%d, j=%d, error=%f.\n", i, j, error);
            validity=0;
        end
    end
end

%% Frobenius Reciprocity
for i=1:s
    for j=1:s
        for k=1:s
            d = M{i}(k,j); % N_{ij}^k
            if M{j}(dual(1,i),dual(1,k)) ~= d 
                fprintf("N_{jk^*}^{i^*} ~= N_{ij}^k; i=%d, j=%d, k=%d, i^*=%d, j^*=%d, k^*=%d.\n", i,j,k,dual(1,i),dual(1,j),dual(1,k));
                validity=0;
            end
            if M{dual(1,k)}(dual(1,j),i) ~= d 
                fprintf("N_{k^* i}^{j^*} ~= N_{ij}^k; i=%d, j=%d, k=%d, i^*=%d, j^*=%d, k^*=%d.\n", i,j,k,dual(1,i),dual(1,j),dual(1,k));
                validity=0;
            end
            if M{k}(i,dual(1,j)) ~= d 
                fprintf("N_{kj^*}^{i} ~= N_{ij}^k; i=%d, j=%d, k=%d, i^*=%d, j^*=%d, k^*=%d.\n", i,j,k,dual(1,i),dual(1,j),dual(1,k));
                validity=0;
            end
            if M{dual(1,i)}(j,k) ~= d 
                fprintf("N_{i^*k}^{j} ~= N_{ij}^k; i=%d, j=%d, k=%d, i^*=%d, j^*=%d, k^*=%d.\n", i,j,k,dual(1,i),dual(1,j),dual(1,k));
                validity=0;
            end
            if M{dual(1,j)}(dual(1,k), dual(1,i)) ~= d 
                fprintf("N_{j^* i^*}^{k^*} ~= N_{ij}^k; i=%d, j=%d, k=%d, i^*=%d, j^*=%d, k^*=%d.\n", i,j,k,dual(1,i),dual(1,j),dual(1,k));
                validity=0;
            end
        end
    end
end
end
%% Exclusion of fusion rings using the NC Minkowski inequality
for multiplicity=1:16

%loads Fusion Ring Dataset 'DSL' with given multiplicity 
%DSL{i}{j}: jth fusion ring with rank i+1
FusionRingDataSet_load_DSL
[a,len_DSL]=size(DSL); %a=1

%Initialization
% Comment the first two lines and uncomment the last two,
% only if running the code for the first time.
load('FusionRingDataSet_Solutions.mat','FusionRingDataSet_Solutions')
load('FusionRingDataSet_Summary.mat','FusionRingDataSet_Summary')
%FusionRingDataSet_Solutions{multiplicity}=cell(1,len_DSL);
%FusionRingDataSet_Summary{multiplicity}=zeros(4,len_DSL);



for i=1:len_DSL %i=rank-1
tic;

[a,len_DSL_i]=size(DSL{i}); %a=1, len_DSL_i=# rank i+1 fusion rings
fprintf("The number of classes of fusion rings in DSL: %d \n", len_DSL);
fprintf("The number of fusion rings in the %dth class: %d \n", i, len_DSL_i);

%% Checking the multiplicity (optional)
%{
max_coef=0;
for j=1:len_DSL_i % DSL{i}{j} = fusion ring
    [a,len_DSL_ij]=size(DSL{i}{j}); %a=1
    max_coef_ij =0;
    for k=1:len_DSL_ij % DSL{i}{j}{k} = fusion matrix
        m=max(max(DSL{i}{j}{k}));
        if m > max_coef_ij
            max_coef_ij=m;
        end
    end
    if max_coef_ij > max_coef
        max_coef = max_coef_ij;
    end
end

fprintf("Maximum fusion coefficient: %d\n", max_coef)
%}



%% Checking the validity of the fusion rings (optional)
%{
for j=1:len_DSL_i % DSL{i}{j} = fusion ring
    if FR_check(DSL{i}{j})
        %fprintf("DSL{%d}{%d} is a valid fusion rule.\n", i, j)
    else
        fprintf("ERROR: DSL{%d}{%d} is not a valid fusion rule.\n", i, j)
    end
end
%}

%% Minkowski test
% Implements the gradient descent method to obtain small RHS/LHS
% where RHS, LHS are the right- and left-hand side of the Minkowski ineq.

%i = rank - 1; %rank = i+1, 
%[a,len_DSL_i]=size(DSL{i}); %a=1

iterations = 20; % # of times to apply gradient descent
nonunitary=zeros(1,len_DSL_i); % will record the index of excluded fusion rings here.
FusionRingDataSet_Solutions{multiplicity}{i}=cell(1,len_DSL_i); % will record inputs violating ineq. here.
for j=1:len_DSL_i % DSL{i}{j} = jth fusion ring with rank i+1
    fprintf("DSL{%d}{%d}\n", i, j)
    LessThan1_i_j=cell(0,0); % will record inputs violating ineq. here.
    count_sol = 0; % counts inputs that fails Minkowski ineq.
    for n = 1:iterations
        %DSL{i}{j}
        [is_unitary,V] = FR_Euler_method( DSL{i}{j}, 0.0001, 0.05, 1000 ); 
        % FR_Euler_method( fusion_ring , Delta_x, Delta_t, steps )
        if is_unitary == 0
            nonunitary(1,j)=1;
            LessThan1_i_j=[LessThan1_i_j;V];
            [a,b]=size(V);
            count_sol = count_sol+a;
            if count_sol > 20
                break
            end
        end
    end
    % Saving the vectors that fail to pass Minkowski Criterion.
    FusionRingDataSet_Solutions{multiplicity}{i}{1,j}=LessThan1_i_j;
    save('FusionRingDataSet_Solutions.mat','FusionRingDataSet_Solutions')
end

excluded=sum(nonunitary);
fprintf("Excluded fusion rings: %d out of %d \n", excluded,len_DSL_i);
for j= 1:len_DSL_i
    if nonunitary(1,j)==1
        fprintf("%d   ",j);
    end
end
fprintf("\n\n")






FusionRingDataSet_Summary{multiplicity}(1,i)=len_DSL_i; % number of fusion rings
FusionRingDataSet_Summary{multiplicity}(2,i)=excluded; % number of excluded fusion rings
FusionRingDataSet_Summary{multiplicity}(3,i)=excluded/len_DSL_i; % exclusion rate
FusionRingDataSet_Summary{multiplicity}(4,i)=toc; % running time
save('FusionRingDataSet_Summary.mat','FusionRingDataSet_Summary')

FusionRingDataSet_Solutions
FusionRingDataSet_Summary{multiplicity}



end
end
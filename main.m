% OCT. 4. 2016. YEJIN KIM
% main.m : load data and run the phenotype code



clear;

% tensor_toolbox should be downloaded in advance
addpath(genpath('./tensor_toolbox')); 

% the number of candidate phenotypes
rank=50;

%% LOAD COUNT DATA
fileName= 'count.csv'; % count, subject_id, diagnosis(icd9_3), prescription
count_mat=csvread(fileName); 
sz = max(count_mat);


%% LOAD SIMILARITY MATRIX
isSimMatAvail = 0;
if (isSimMatAvail)
    D = dlmread('similarities.txt', ' '); %pairwise similarity matrix
    D(D < 0) = 0; % ignore dissimilarity score if exists
    
    only consider k nearest similar items
    k_near = floor(log2(sz(3))) +1;
    for i=1:sz(3)
        [~, idx]=sort(D(i,1:sz(3)), 'descend');
        D(i, setdiff(1:sz(3), idx(1:(1+k_near)))) = 0;
        
    end
    
    for i=(sz(3)+1):(sz(3) + sz(4))
        [~, idx]=sort(D(i, (sz(3)+1):(sz(3) + sz(4))), 'descend');
        idx = idx + sz(3);
        D(i, setdiff((sz(3)+1):(sz(3) + sz(4)), idx(1:(1+k_near)))) = 0;
        
    end
    
    D = (D + D')/2; % symmetirc
    
    diagnoal_mat = diag(sum(D,2));
    D=((sqrt(diagnoal_mat))\ D) /(sqrt(diagnoal_mat)); %Normalized cut similarity
    
else
    % orthogonality
    D = eye(size(sz(3)+sz(4)));
end

%% MAKE OBSERVED TENSOR
O=sptensor(count_mat(:, 2:4), count_mat(:, 1), sz(2:end)); %index, count, dimension

%% Training set and label
% load label
dead = csvread('label.csv'); % subject_id (cd), label (1 if death, -1 o/w)


label=-ones(sz(2),1);
label(dead) = 1;

omega = 0;
mu = 0;

% training set
indices = crossvalind('Kfold', sz(2), 10);
%parfor cv = 1:10
for cv = 1:1
    test = (indices == cv); 
    train = ~test;
    
    runPhenotyping(O, D, label, rank, mu, omega, train, test, cv)
end






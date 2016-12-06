%  OCT 4 2016 VERSION14 YEJIN KIM
% INPUT
% O: observed 3-way tensor
% D: similarity symmetric matrix of size (numMed + numLab) by (numMed + numLab)
% label: outcomes
% rank: rank of CP factorization
% mu: 


% OUTPUT
% T: ktensor
 

function runPhenotyping(O, D, label, rank, mu, omega, train, test, cv)

dim = size(O);
maxiter = 100;
N=3;
mumax=1e10;
omegamax= 1e10;

update_ratio=1.2;

sparsity=[0.000001, 0.001, 0.001];

normO = norm(O);
fitchangetol = 0.0001;

%% initialization

% B
B = cell(1,N);

if (mu ~= 0 && ~(abs(sum(D(1, 1:dim(2))) -1) < 1e-6)) %initial B is from spectral clustering   
    B{1} = rand(dim(1), rank);
    B{1}(B{1} <= sparsity(1)) = 0;
    
    Dnn = D(1:dim(2), 1:dim(2));
    B{2} = spectralClustering(Dnn, dim(2), rank);
    
    Dnn = D((dim(2)+1):(dim(2) + dim(3)), (dim(2)+1):(dim(2) + dim(3)));
    B{3} = spectralClustering(Dnn, dim(3), rank);
    
    fileName= strcat('omega', num2str(omega), 'mu', num2str(mu), 'cv', num2str(cv),  '.mat');
    
else %if no similarity information is available (e.g. only with orthogonality matrix)
    for n = 1:N
        B{n} = rand(dim(n), rank);
        B{n}(B{n} <= sparsity(n)) = 0;
        
    end
    
    fileName= strcat('omega', num2str(omega), 'mu', num2str(mu), 'noInitial', 'cv', num2str(cv),  '.mat');
end


[Coef, FitInfo] = lassoglm(B{1}, label > 0, 'binomial', 'CV', 5, 'Alpha', 0.1);
Theta = Coef(:, FitInfo.IndexMinDeviance);
theta = FitInfo.Intercept(FitInfo.IndexMinDeviance);


%% main loop
fit = 0;
for iter = 1 : maxiter
      
    fitold = fit;
    
    for n = 1:N
        % update B{n}

        
        
        pitpi=ones(rank, rank);
        for i = [1:n-1, n+1:N]
            pitpi = pitpi .* (B{i}' * B{i});
        end       

        term1 = 2 * B{n} * pitpi - 2 * mttkrp(O, B, n);
        term1_H = 2 * kron(pitpi, speye(dim(n)));
        
        if n ~= 1
            
            if (mu ~= 0)
                d_idx1= sum(dim(2:(n-1)))+1:sum(dim(2:n));
                d_idx2= d_idx1;

                Dnn = D(d_idx1, d_idx2);

                term2 = 4 * mu * (B{n} * B{n}' - Dnn) * B{n};
                
         
                term2_H = 4 * mu * (2 * B{n}(:) * B{n}(:)' + (B{n}(:)' * B{n}(:))*speye(dim(n)*rank) - kron(speye(rank), Dnn) );
                
                term3 = 0;

            else
                term2 = 0; term3 = 0;
                term2_H = 0;
            end

            
            
            gradient = term1 + term2 + term3; 
            H = term1_H + term2_H;
  
            
            Bnew = B{n}(:) - H\gradient(:);
            B{n} = reshape(Bnew, [dim(n), rank]);
            
            idx = B{n} < 0;
            if (any(idx(:) == 1))
                B{n}(idx) = 0;
            end
            
        else

            % discriminative

            if omega~=0
                for i = 1:dim(1)
                    if train(i) == 1
                        
                        
                        gradient_row = term1(i, :) - omega * label(i)/(dim(1)*0.9)  *  1/(1+ exp(label(i)*(B{n}(i, :) * Theta + theta))) * Theta';                       
                        
                        
                        H_row = 2 * pitpi + omega /(dim(1)*0.9) * 1/ (2 + exp((B{n}(i, :) * Theta + theta)) + exp(-(B{n}(i, :) * Theta + theta)) ) * (Theta * Theta');

                        
                        B{n}(i, :) = B{n}(i, :) - (H_row \ gradient_row')';

                    else
                        
                        gradient_row =  term1(i, :);
                        H_row = 2 * pitpi;
                        

                        B{n}(i, :) = B{n}(i, :) - (H_row \ gradient_row')';
                        
                    end
                end
                
                idx = B{n} < 0;
                if (any(idx(:) == 1))
                    B{n}(idx) = 0;
                end
                
                
                % update Theta, theta
                [Coef, FitInfo] = lassoglm(B{1}, label > 0, 'binomial', 'CV', 5, 'Alpha', 0.1);
                Theta = Coef(:, FitInfo.IndexMinDeviance);
                theta = FitInfo.Intercept(FitInfo.IndexMinDeviance);

                
            else % if omega = 0
                gradient = term1;
                
                H = 2 * kron(pitpi, speye(dim(n)));
                
                Bnew = B{n}(:) - H\gradient(:);
                B{n} = reshape(Bnew, [dim(n), rank]);
                
                idx = B{n} < 0;
                if (any(idx(:) == 1))
                    B{n}(idx) = 0;
                end
                
            end
            
            
        end

    end
    
    %% update parameter
    mu = min(update_ratio * mu, mumax);
    omega = min(update_ratio * omega, omegamax);

    
    T = ktensor(B);

    normresidual = sqrt( normO^2 + norm(T)^2 - 2 * innerprod(O,T) );
    
    fit = 1 - (normresidual / normO);
    
    fitchange = abs(fitold - fit);
    
    rmse=normresidual/(sqrt(nnz(O)));
    
    fprintf(' Iter %2d: fit = %e fitdelta = %7.1e RMSE = %e\n', iter, fit, fitchange, rmse);
    
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end
    
    
end




%% clean up
T=arrange(T);


%% Discriminativeness
% all
[Coef, FitInfo] = lassoglm(T.U{1}, label > 0, 'binomial', 'CV', 5, 'Alpha', 0.1);
Theta = Coef(:, FitInfo.IndexMinDeviance);
theta = FitInfo.Intercept(FitInfo.IndexMinDeviance);

ypred = glmval([theta; Theta], T.U{1}(test, :),'probit');
[~, ~, ~, auc, optrocpt]=perfcurve(label(test), ypred, '1');
[~, ~, ~, inmodel ,~, ~, ~] = stepwisefit(T.U{1}(train,:), label(train, :));

% selected features
[Coef, FitInfo] = lassoglm(T.U{1}(:, inmodel), label > 0, 'binomial', 'CV', 5, 'Alpha', 0.1);
Theta_select = Coef(:, FitInfo.IndexMinDeviance);
theta_select = FitInfo.Intercept(FitInfo.IndexMinDeviance);
ypred = glmval([theta_select; Theta_select], T.U{1}(test, inmodel),'probit');
[~, ~, ~, auc_select, optrocpt_select]=perfcurve(label(test), ypred, '1');


%% Sparseness
for i = 1: N
    T{i}(T{i} <= sparsity(i)) = 0; % sparse factor
    
end

numZeros = zeros(1,N);
for k=1:N
    numZeros(k) = sum(T.U{k}(:)==0);
end

sumValues = zeros(1,N);
for k=1:N
    sumValues(k) = sum(T.U{k}(:));
end


total = dim(2)*rank + dim(3)*rank;
len_per_phe=(total - (numZeros(2) + numZeros(3)))/rank;
avg_ratio_per_phe = len_per_phe/(dim(2) + dim(3));

%average overlap
avgOverlap =  (sum(1-pdist(T.U{2}', 'cosine' )) + sum(1-pdist(T.U{3}', 'cosine')) )/ (rank*(rank-1));

save(fileName, 'T', 'B', 'Theta', 'Theta_select', 'train', 'test', 'rmse', 'label', 'auc', 'auc_select', 'optrocpt', 'optrocpt_select', 'inmodel', 'avg_ratio_per_phe', 'numZeros','avgOverlap');



end

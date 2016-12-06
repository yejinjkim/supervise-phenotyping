function B = spectralClustering(Dnn, sz, rank)

%Unnormalized graph laplacian
L = diag(sum(Dnn,2)) - Dnn;
diagnoal_mat = diag(sum(Dnn,2));

[eVec, ~]=eig(L, diagnoal_mat);

U = eVec(:, 1:rank);

idx = kmeans(U, rank);

B = zeros(sz, rank);
for r = 1:rank
    B(:,r) = (idx ==r);
end

end

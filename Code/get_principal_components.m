function pc = get_principal_components(data,Zscore)
% Get principal components (see Stock & Watson 2002, eq 6)
% data is (T, nvar)

[T,n] = size(data);

if Zscore == 1
      for i=1:n
            data(:,i) = zscore(data(:,i));
      end
end

% Get VC matrix
sigma = 1/T*data'*data;

[V,lamb] = eig(sigma); % diag(lamb) are the eigenvalues, V the eigenvectors
V = real(V);
lamb = real(lamb);
pc = data*V./n;

end

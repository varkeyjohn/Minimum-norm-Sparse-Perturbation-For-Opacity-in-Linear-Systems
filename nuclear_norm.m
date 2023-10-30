function obj = nuclear_norm(X,n)
    if size(X, 1) > 2000
        sigma = svds(X, 100);
    else
        sigma = svds(X,n);
    end
    obj = sum(sigma);
end
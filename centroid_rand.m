function [ c ] = centroid_rand( stride, supp, w )
% generate random centroid using mvnrnd()
  n = length(stride);
  m = length(w);
    
  avg_stride = ceil(mean(stride));  

  % Compute random initial guess
    c_means = supp * w' / n;
    zm = supp - repmat(c_means, [1, m]);
    c_covs = zm * diag(w) * zm' / n;
    c.supp = mvnrnd(c_means', c_covs, avg_stride)';
    %c.w = rand(1,avg_stride); c.w = c.w/sum(c.w);
    c.w = 1/avg_stride * ones(1,avg_stride);


end


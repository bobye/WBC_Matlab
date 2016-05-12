function [ c, resampler ] = centroid_rand( stride, supp, w, options )
% generate random centroid using mvnrnd() [default] or kmeans()
  n = length(stride);
  m = length(w);
  dim = size(supp,1);

  if isfield(options, 'support_size')
      support_size=options.support_size;
  else
      support_size = ceil(mean(stride));  
  end
  
  init_method='kmeans';
  if isfield(options, 'init_method');
      init_method=options.init_method;
  end

  if strcmp(init_method, 'mvnrnd')
  % Compute random initial guess
    c_means = supp * w' / (n-dim);
    zm = supp - repmat(c_means, [1, m]);
    c_covs = zm * diag(w) * zm' / (n-dim) + 1e-4 * eye(dim);
    c.supp = mvnrnd(c_means', c_covs, support_size)';
    %c.w = rand(1,avg_stride); c.w = c.w/sum(c.w);
    c.w = 1/support_size * ones(1,support_size);
  elseif strcmp(init_method, 'kmeans')
    [idx, c.supp] = kmeans(supp', support_size);
    c.supp = c.supp';
    c.w = histc(idx, 1:support_size);
    c.w = c.w' / sum(c.w);
  end
    function P = onesupp (k)
        P = mvnrnd(c_means',c_covs, k)';
    end

   if strcmp(init_method, 'mvnrnd')
       resampler = @(k) onesupp(k);
   else
       resampler = [];
   end
end


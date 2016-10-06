function [c] = centroid_sphEM_GMM(stride, supp, w, c0, options) 
  if isempty(c0)
    error('Please give a GMM barycenter initialization.');
  else
    c=c0;
  end
  support_size=length(c.w);
  d = floor(sqrt(size(supp,1)));
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  X = zeros(support_size, m);
  D = zeros(n,1);  

  % create buffering data
  XX = cell(n,1);
  suppx = cell(n,1);
  wx = cell(n,1);
  strips=cell(n,1);
  for iter=1:n
      strips{iter} = posvec(iter):(posvec(iter)+stride(iter)-1);
      suppx{iter} = supp(:,strips{iter});
      wx{iter} = w(strips{iter});
  end
  
  nIter = 100;     
  if isfield(options, 'max_iters')
      nIter=options.badmm_max_iters;
  end
  
  fval = Inf;
  for iter = 1:nIter
      fval0 = fval;
      C = pdist2(c.supp(1:d,:)', supp(1:d,:)', 'sqeuclidean');
      C = C + gaussian_wd(c.supp((d+1):end,:), supp((d+1):end,:));
      for it=1:n
        [D(it), XX{it}, ~] = kantorovich_with_cost(C(:, strips{it}), c.w, wx{it}); 
      end
      fval=mean(D);
      for j=1:n
        X(:,strips{j}) = XX{j};
      end
      c.supp(1:d,:) = supp(1:d,:) * X' ./ repmat(sum(X,2)', [d, 1]);
      c.supp((d+1):end,:) = gaussian_mean(supp((d+1):end,:), X, c.supp((d+1):end,:));
      fprintf('\t %d %f\n', iter, fval);
      if abs(fval - fval0) < 1E-6 * fval
        break;
      end

  end
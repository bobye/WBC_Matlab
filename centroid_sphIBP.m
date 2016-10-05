function [c] = centroid_sphIBP(stride, supp, w, c0, options) 
% Single phase centroid using Iterative Bregman Projection [Benamou et. al. 2015]
% This code has been created by Jianbo Ye (jxy198 [AT] ist.psu.edu).
  d = size(supp,1);
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];

  if isempty(c0)
    c=centroid_init(stride, supp, w, options);
  else
    c=c0;
  end
  support_size=length(c.w);

  %save cstart.mat
  %load(['cstart' num2str(n) '-' num2str(avg_stride) '.mat']);
  %return;
  

  spIDX_rows = zeros(support_size * m,1);
  spIDX_cols = zeros(support_size * m,1);
  for i=1:n
      [xx, yy] = meshgrid((i-1)*support_size + (1:support_size), posvec(i):posvec(i+1)-1);
      ii = support_size*(posvec(i)-1) + (1:(support_size*stride(i)));
      spIDX_rows(ii) = xx';
      spIDX_cols(ii) = yy';
  end
  spIDX = repmat(speye(support_size), [1, n]);
  
  % initialization
  C = pdist2(c.supp', supp', 'sqeuclidean');
  
  nIter = 20000;
  if isfield(options, 'ibp_max_iters')
      nIter = options.ibp_max_iters;
  end
  
  if isfield(options, 'ibp_vareps')
      rho = options.ibp_vareps * median(median(pdist2(c.supp', supp', 'sqeuclidean')));
  else
      rho = 0.01 * median(median(pdist2(c.supp', supp', 'sqeuclidean')));
  end
  
  if isfield(options, 'ibp_tol')
      ibp_tol = options.ibp_tol;
  else
      ibp_tol = 1E-4; % no updates of support
  end
  
  
  xi=exp(-C / rho);
  xi(xi<1e-200)=1e-200; % add trick to avoid program breaking down
  xi=sparse(spIDX_rows, spIDX_cols, xi(:), support_size * n, m);
  v = ones(m, 1);
  w1=w';
  fprintf('\n');
  obj=Inf;
  tol=Inf;
  for iter = 1:nIter
    w0=repmat(c.w', n, 1);
    u=w0 ./ full(xi*v);
    v=w1 ./ full(xi'*u);
    c.w = geomean(reshape(u .* full(xi * v), support_size, n), 2)';
    
    if (mod(iter, 10) == 0)
    tol = norm(reshape(full(spdiags(u, 0, n*support_size, n*support_size) * xi * spdiags(v, 0, m, m) * ones(m,1)) ...
        - w0, support_size, n), Inf); 
    end
    

    if tol < ibp_tol && ~isfield(options, 'support_points')
        c_back = c;
        X=full(spIDX * spdiags(u, 0, support_size*n, support_size*n) * xi * spdiags(v, 0, m, m));
        c.supp = supp * X' ./ repmat(sum(X,2)', [d, 1]);
        C = pdist2(c.supp', supp', 'sqeuclidean');
        xi=exp(-C / rho);
        xi(xi<1e-200)=1e-200; % add trick to avoid program breaking down
        xi=sparse(spIDX_rows, spIDX_cols, xi(:), support_size * n, m);
        v = ones(m, 1);
        last_obj=obj;
        obj=sum(C(:).*X(:))/n;
        fprintf('\t %d %f\n', iter, obj);   
        if (obj>last_obj)
            c = c_back;
            fprintf('terminate!\n');
            break;
        end
        tol=Inf;
    end
    
    if (tol < ibp_tol && isfield(options, 'support_points'))
        fprintf('iter = %d\n', iter);
        break;
    end

  end
  
end
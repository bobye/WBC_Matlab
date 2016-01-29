function [c] = centroid_sphIBP(stride, supp, w, c0) 
  d = size(supp,1);
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  avg_stride = ceil(mean(stride));  

  if isempty(c0) || length(c0.w)~=avg_stride
    c=centroid_rand(stride, supp, w);
  else
    c=c0;
  end
  
  %centroid_sphEnergy(stride, supp, w, c);
  
  %save cstart.mat
  %save(['cstart' num2str(n) '.mat'], 'c', 'avg_stride');
  load(['cstart' num2str(n) '.mat']);
  %return;
  

  spIDX_rows = zeros(avg_stride * m,1);
  spIDX_cols = zeros(avg_stride * m,1);
  for i=1:n
      [xx, yy] = meshgrid((i-1)*avg_stride + (1:avg_stride), posvec(i):posvec(i+1)-1);
      ii = avg_stride*(posvec(i)-1) + (1:(avg_stride*stride(i)));
      spIDX_rows(ii) = xx';
      spIDX_cols(ii) = yy';
  end
  spIDX = repmat(speye(avg_stride), [1, n]);
  
  % initialization
  C = pdist2(c.supp', supp', 'sqeuclidean');
  
  nIter = 2000;
  rho = 1E-2 * mean(mean(pdist2(c.supp', supp', 'sqeuclidean')));
  
  xi=exp(-C / rho);
  xi=sparse(spIDX_rows, spIDX_cols, xi(:), avg_stride * n, m);

  v = ones(m, 1);
  w1=w';
  for iter = 1:nIter
    w0=repmat(c.w', n, 1);
    u=w0 ./ full(xi*v);
    v=w1 ./ full(xi'*u);
    %c.w = mean(reshape(u .* full(xi * v), avg_stride, n), 2)';
    c.w = geomean(reshape(u .* full(xi * v), avg_stride, n), 2)';
    if false && mod(iter, 10)==0
        X=full(spIDX * spdiags(u, 0, avg_stride*n, avg_stride*n) * xi * spdiags(v, 0, m, m));
        c.supp = supp * X' ./ repmat(sum(X,2)', [d, 1]);
        C = pdist2(c.supp', supp', 'sqeuclidean');
        xi=exp(-C / rho);
        xi=sparse(spIDX_rows, spIDX_cols, xi(:), avg_stride * n, m);
    end
    
    % output
    if false && mod(iter, 100) == 0
        fprintf('\t %d %f\n', iter, sum(C(:).*X(:))/n);        
    end    
  end
  
end
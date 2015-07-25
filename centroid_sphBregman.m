function [c] = centroid_sphBregman(stride, supp, w, c0) 
  % The algorithmic prototype of Wasserstein Barycenter using Bregman ADMM
  % This approach has been described in the following paper:
  %     Jianbo Ye, Panruo Wu, James Z. Wang and Jia Li, Fast Discrete
  %     Distribution Clustering Under Wasserstein Distance, submitted 2015
  % 
  % This code has been created by Jianbo Ye (jxy198 [AT] ist.psu.edu).
  % 
  dim = size(supp,1);
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  avg_stride = ceil(mean(stride));  

  if isempty(c0) || length(c0.w)~=avg_stride
    c=centroid_rand(stride, supp, w);
  else
    c=c0;
  end
  
  %save cstart.mat
  save(['cstart' num2str(n) '.mat'], 'c', 'avg_stride');
  %return;
  X = zeros(avg_stride, m);
  Y = zeros(size(X)); Z = X;
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
  for i=1:n
      Z(:,posvec(i):posvec(i+1)-1) = 1/(avg_stride*stride(i));
  end
  C = pdist2(c.supp', supp', 'sqeuclidean');
  
  nIter = 2000;     
  rho = 2.*mean(mean(pdist2(c.supp', supp', 'sqeuclidean')));
  for iter = 1:nIter
      % update X
      X = Z .* exp((C+Y)/(-rho)) + eps;
      X = bsxfun(@times, X', w'./sum(X)')';
      
      % update Z
      Z0 = Z;
      Z = X .* exp(Y/rho) + eps;
      spZ = sparse(spIDX_rows, spIDX_cols, Z(:), avg_stride * n, m);
      tmp = full(sum(spZ, 2)); tmp = reshape(tmp, [avg_stride, n]);
      dg = bsxfun(@times, 1./tmp, c.w'); 
      dg = sparse(1:avg_stride*n, 1:avg_stride*n, dg(:));
      Z = full(spIDX * dg * spZ);
      
      % update Y      
      Y = Y + rho * (X - Z);
      
      % update c.w
      tmp = bsxfun(@times, tmp', 1./sum(tmp)');
      sumW = sum(tmp);
      c.w = sumW / sum(sumW);
      
      % update c.supp and compute C (lazy)
      if mod(iter, 10)==0
        c.supp = supp * X' ./ repmat(sum(X,2)', [dim, 1]);      
        C = pdist2(c.supp', supp', 'sqeuclidean');
      end
      
      % The constraint X=Z are not necessarily strongly enforced
      % during the update of w, which makes it suitable to reset
      % lagrangian multipler after a few iterations
      if (mod(iter, 100) == 0)
%          Y(:,:) = 0;
%           if primres > 10*dualres
%             rho = 2 * rho;
%             fprintf(' *2');
%           elseif 10*primres < dualres
%             rho = rho / 2;
%             fprintf(' /2');
%           end
      end
      
      % output
      if (mod(iter, 100) == 0)
          primres = norm(X-Z,'fro')/norm(Z,'fro');
          dualres = norm(Z-Z0,'fro')/norm(Z,'fro');
          fprintf('\t %d %f %f %f ', iter, sum(sum(C.*X))/n, ...
              primres, dualres);
          fprintf('\n');          
      end
  end
   
end
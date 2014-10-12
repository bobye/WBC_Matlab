function [c] = centroid_sphBregman(stride, supp, w, c0) 

  dim = size(supp,1);
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  posStride = cell(n,1); 
  for i=1:n
      posStride{i} = posvec(i):posvec(i+1)-1;
  end
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
  Y = X; Z = X; C = X;
  %D = zeros(n,1);
  
  
  % initialization
  for i=1:n
      Z(:,posStride{i}) = 1/(avg_stride*stride(i));
  end
  C = pdist2(c.supp', supp', 'sqeuclidean');
  
  nIter = 2000;     
  rho = mean(mean(pdist2(c.supp', supp', 'sqeuclidean')));
  for iter = 1:nIter
      % update X
      X = Z .* exp(- (C+Y)/rho);
      X = X .* repmat(w./sum(X), [avg_stride,1]);
      
      % update Z
      sumlogW = zeros(1,avg_stride);
      Z0 = Z;
      Z = X .* exp(Y/rho);
      sumW = sum(Z,2)';
      for i=1:n
          tmp = sum(Z(:,posStride{i}),2)';
          sumlogW = sumlogW + log(tmp).* tmp;
          Z(:,posStride{i}) = bsxfun(@times, Z(:,posStride{i})', c.w./tmp)'; % MATLAB is slow
      end
      % update Y      
      Y = Y + rho * (X - Z);      
      % update c.w
      sumW = exp(sumlogW ./ sumW);
      c.w = sumW / sum(sumW);
      % update c.supp and compute C (lazy)
      if mod(iter-1, 100)==0
        c.supp = supp * X' ./ repmat(n*c.w, [dim, 1]);      
        C = pdist2(c.supp', supp', 'sqeuclidean');
      end
      
      % The constraint X=Z are not necessarily strongly enforced
      % during the update of w, which makes it suitable to reset
      % lagrangian multipler after a few iterations
      if (mod(iter, 10) == 0)
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
          primres = norm(X(:)-Z(:))/norm(Z(:));
          dualres = norm(Z(:)-Z0(:))/norm(Z(:));
          fprintf('\t %d %f %f %f ', iter, sum(sum(C.*X))/n, ...
              primres, dualres);
          fprintf('\n');          
      end
  end
   
end
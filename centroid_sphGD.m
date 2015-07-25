function [ c ] = centroid_sphGD( stride, supp, w, c0)
  % The algorithmic prototype of Wasserstein Barycenter using subgradient
  % descent method with re-parametrization.
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
  
  % compute obj and grad
  function  [obj grad] = d2energy(c_supp, c_w)
  
    if (any(c_w < 0)) 
      c_w(c_w>0) = 0; c_w = c_w / sum(c_w);
      %fprintf('%f ',c_w);
      %error ('w can not be smaller than zero');  
    end
  
    grad = zeros(avg_stride,1);
    for it=1:n                              
      [D(it), XX{it}, lambda] = kantorovich(c_supp, c_w, suppx{it}, wx{it});      

      f=lambda.eqlin(1:avg_stride); f=f-sum(f)/avg_stride;
      grad = grad - f;
    end
    obj = mean(D);
    grad = grad/n;
    %fprintf(stdoutput, '\n%e', obj );  
  end
  
  nIter = 20;
  fval = Inf;
  for iter = 1:nIter
    tic;
    fval0 = fval;
    [fval, gw] = d2energy(c.supp, c.w);
    fprintf('\t%e ', fval ); 
    for j=1:n
        X(:,strips{j}) = XX{j};
    end
    c.supp = supp * X' ./ repmat(n*c.w, [dim, 1]);
    c.w = c.w .* exp(-8E-3 * gw' .* c.w .* (1-c.w));
    c.w = c.w / sum(c.w);
    toc;
    if abs(fval - fval0) < 1E-6 * fval
        break;
    end    
  end
end


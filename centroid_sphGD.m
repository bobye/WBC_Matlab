function [ c ] = centroid_sphGD( stride, supp, w, c0, options)
  % The algorithmic prototype of Wasserstein Barycenter using subgradient
  % descent method with re-parametrization.
  % 
  % This code has been created by Jianbo Ye (jxy198 [AT] ist.psu.edu).
  %   
  dim = size(supp,1);
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
  %save(['cstart' num2str(n) '.mat'], 'c', 'support_size');
  %return;

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
  
  % compute obj and grad
  function  [obj grad] = d2energy(c_supp, c_w)
  
    if (any(c_w < 0)) 
      c_w(c_w>0) = 0; c_w = c_w / sum(c_w);
      %fprintf('%f ',c_w);
      %error ('w can not be smaller than zero');  
    end
  
    grad = zeros(support_size,1);
    for it=1:n                              
      [D(it), XX{it}, lambda] = kantorovich(c_supp, c_w, suppx{it}, wx{it});      

      f=lambda.eqlin(1:support_size); f=f-sum(f)/support_size;
      grad = grad - f;
    end
    obj = mean(D);
    grad = grad/n;
    %fprintf(stdoutput, '\n%e', obj );  
  end
  
  nIter = 20;
  if isfield(options, 'max_iters')
      nIter = options.max_iters;
  end
  step_size=0.008;
  if isfield(options, 'gd_step_size')
      step_size=options.gd_step_size;
  end
  
  fval = Inf;
  for iter = 1:nIter
    fval0 = fval;
    [fval, gw] = d2energy(c.supp, c.w);
    fprintf('\n%d\t%e ', iter, fval ); 
    for j=1:n
        X(:,strips{j}) = XX{j};
    end
    if ~isfield(options, 'support_points')
        c.supp = supp * X' ./ repmat(n*c.w, [dim, 1]);
    end
    c.w = c.w .* exp(-step_size * gw' .* c.w .* (1-c.w));
    c.w = c.w / sum(c.w);
    if abs(fval - fval0) < 1E-6 * fval
        break;
    end    
  end
end


function [c] = centroid_sphADMM(stride, supp, w, c0, options) 
  % Single phase centroid using ADMM
  % The algorithmic prototype of Wasserstein Barycenter using ADMM
  % This approach has been described in the following paper:
  %     Jianbo Ye and Jia Li, Scaling Up Discrete 
  %     Distribution Clustering Using ADMM, ICIP 2014
  % 
  % This code has been created by Jianbo Ye (jxy198 [AT] ist.psu.edu).
  
if isfield(options, 'mosek_path')
    addpath(options.mosek_path);
end
  
% Re-prepare
  global A B;
  global stdoutput qpoptim_options;
  
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

  %save(['cstart' num2str(n) '.mat'], 'c', 'avg_stride');
  
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
  
  

  function  obj = d2energy(warm)
  for it=1:n                              
    if warm
      [D(it), XX{it}] = kantorovich(c.supp, c.w, suppx{it}, wx{it}, XX{it});
    else
      [D(it), XX{it}] = kantorovich(c.supp, c.w, suppx{it}, wx{it});        
    end
  end
    obj = mean(D);
    fprintf(stdoutput, '\n\t\t %d\t %f', iter, obj );  
  end

  d2energy(false);


% ADMM optimization

  nIter = 50; 
  admmIter = 10;
  if isfield(options, 'admm_max_iters') && ~isfield(options, 'support_points')
      nIter = options.admm_max_iters;
  elseif isfield(options, 'support_points')
      nIter = 1;
      admmIter=50;
  end 
  suppIter = 1;
  
  if isfield(options, 'admm_inner_iters')
      admmIter=options.admm_inner_iters;
  end
  
  rho0=50;
  if isfield(options, 'admm_rho')
      rho0 = options.admm_rho;
  end

  fprintf(stdoutput,'\n');  
  statusIter = zeros(nIter,1);
  elapsedTime = zeros(nIter,1);
  iter_tol = 1E-6;
  tic; 
  for iter=1:nIter    
    if ~isfield(options, 'support_points')  
    for xsupp=1:suppIter
      % update c.supp
      for j=1:n
        X(:,strips{j}) = XX{j};
      end
      c.supp = supp * X' ./ repmat(sum(X,2)', [dim, 1]);

      % if some components of c.w is zero,
      % we have to reset corresponding components c.supp
      % c.supp(:, abs(c.w)<1E-6) = resampler(sum(abs(c.w)<1E-6));
      
      % setup initial guess for X in ADMM
      % d2energy(true);
    end
    end
    
    % update c.w as well as X, using ADMM

    % empirical parameters
    rho = rho0*mean(D);

    % precompute linear parameters
    C = pdist2(c.supp', supp', 'sqeuclidean');
    Cx = cell(n,1);
    for i=1:n
        Cx{i} = C(:,strips{i});
    end
    
    % lagrange multiplier
    lambda =  zeros(support_size, n); 
    
    for admm=1:admmIter
      % step 1, update X
      
      
      parfor i=1:n
	  vecsize = [support_size * stride(i), 1];
	  
	  x0 = reshape(XX{i}, vecsize);
	  H = rho * B{support_size, stride(i)}; 
	  q = reshape(rho * repmat(lambda(:,i) - c.w', [1, stride(i)]) + Cx{i}, vecsize);
	  Aeq = A{support_size,stride(i)}(support_size+1:end, :);
	  beq = wx{i}';
	  [xtmp] = ... 
	  quadprog(H, q, [], [], Aeq, beq, zeros(vecsize), [], x0, qpoptim_options);
	  XX{i} = reshape(xtmp,[support_size, stride(i)]);
           
      end
      
      for j=1:n
        X(:,strips{j}) = XX{j};
      end
      

      % step 2, update c.w
      w2 = c.w;
      
      H = n*eye(support_size);
      q = - (sum(X, 2) + sum(lambda, 2));
      [c.w] = quadprog(H, q, [], [], ones(1,support_size), 1, zeros(support_size,1), [], c.w', qpoptim_options)';
      
      %H = n * eye(avg_stride) + rho * ones(avg_stride);
      %q = - (sum(X, 2) + sum(lambda, 2) + rho*(1 - mu));
      %[c.w] = quadprog(H, q, [], [], [], [], zeros(avg_stride,1), [], c.w', qpoptim_options)';

      % step 3, update dual variables: lambda and mu
      lambda2 = lambda; 
      
      parfor i=1:n
        lambda(:, i) = lambda(:, i) + sum(XX{i},2) - c.w';
      end
      
      dualres = norm(w2 - c.w);
      primres1 = norm(lambda2 - lambda, 'fro')/sqrt(n*support_size);
      %fprintf(stdoutput, '\t%f\t%f', primres1, dualres);
      
%       if primres1 > 10*dualres
%           pho = 2 * pho;
%           lambda = lambda/2;
%           mu = mu/2;
%           fprintf(stdoutput,' *2');
%       elseif 10*primres1 < dualres
%           pho = pho / 2;
%           lambda = lambda*2;
%           mu = mu*2;
%           fprintf(stdoutput,' /2');
%       end
      % stopping criterion
      if (dualres < 0.005)
          break;
      end
      
    end       

    % sum2one(c.w)
    %c.w = c.w/sum(c.w);

    % output status
    statusIter(iter) = d2energy(false);

    
    elapsedTime(iter) = toc;
    fprintf(stdoutput, '\t%fs', elapsedTime(iter));    
    % pause;
    if iter>1 && abs(statusIter(iter)-statusIter(iter-1))<iter_tol*abs(statusIter(iter))
        break;
    end    
  end
  
  %global statusIterRec;
  %statusIterRec = statusIter;
  
  %h = figure;
  %plot(elapsedTime, statusIter);
  %print(h, '-dpdf', 'centroid_singlephase.pdf');
  
  fprintf(stdoutput, '\n');
  %fprintf(stdoutput, ' %f', c.w);  
  %fprintf(stdoutput, '\n');
  
end


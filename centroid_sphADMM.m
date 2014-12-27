function [c] = centroid_sphADMM(stride, supp, w, c0) 
% Single phase centroid using ADMM
  

  
% Re-prepare
  global A B;
  global stdoutput qpoptim_options;
  
  dim = size(supp,1);
  n = length(stride);
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  avg_stride = ceil(mean(stride));  
  
  
  [c, resampler]=centroid_rand(stride, supp, w);
  if ~isempty(c0) && length(c0.w)==avg_stride
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
  suppIter = 1;
  admmIter = 10;

  fprintf(stdoutput,'\n');  
  statusIter = zeros(nIter,1);
  elapsedTime = zeros(nIter,1);
  tic; 
  for iter=1:nIter    
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
    
    % update c.w as well as X, using ADMM

    % empirical parameters
    pho = 50*mean(D);

    % precompute linear parameters
    C = pdist2(c.supp', supp', 'sqeuclidean');
    Cx = cell(n,1);
    for i=1:n
        Cx{i} = C(:,strips{i});
    end
    
    % lagrange multiplier
    lambda =  zeros(avg_stride, n); 
    
    for admm=1:admmIter
      % step 1, update X
      
      
      parfor i=1:n
	  vecsize = [avg_stride * stride(i), 1];
	  
	  x0 = reshape(XX{i}, vecsize);
	  H = pho * B{avg_stride, stride(i)}; 
	  q = reshape(pho * repmat(lambda(:,i) - c.w', [1, stride(i)]) + Cx{i}, vecsize);
	  Aeq = A{avg_stride,stride(i)}(avg_stride+1:end, :);
	  beq = wx{i}';
	  [xtmp] = ... 
	  quadprog(H, q, [], [], Aeq, beq, zeros(vecsize), [], x0, qpoptim_options);
	  XX{i} = reshape(xtmp,[avg_stride, stride(i)]);
           
      end
      
      for j=1:n
        X(:,strips{j}) = XX{j};
      end
      

      % step 2, update c.w
      w2 = c.w;
      
      H = n*eye(avg_stride);
      q = - (sum(X, 2) + sum(lambda, 2));
      [c.w] = quadprog(H, q, [], [], ones(1,avg_stride), 1, zeros(avg_stride,1), [], c.w', qpoptim_options)';
      
      %H = n * eye(avg_stride) + rho * ones(avg_stride);
      %q = - (sum(X, 2) + sum(lambda, 2) + rho*(1 - mu));
      %[c.w] = quadprog(H, q, [], [], [], [], zeros(avg_stride,1), [], c.w', qpoptim_options)';

      % step 3, update dual variables: lambda and mu
      lambda2 = lambda; 
      
      parfor i=1:n
        lambda(:, i) = lambda(:, i) + sum(XX{i},2) - c.w';
      end
      
      dualres = norm(w2 - c.w);
      primres1 = norm(lambda2 - lambda, 'fro')/sqrt(n*avg_stride);
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
  end
  
  global statusIterRec;
  statusIterRec = statusIter;
  
  %h = figure;
  plot(elapsedTime, statusIter);
  %print(h, '-dpdf', 'centroid_singlephase.pdf');
  
  fprintf(stdoutput, '\n');
  fprintf(stdoutput, ' %f', c.w);  
  fprintf(stdoutput, '\n');
  
end


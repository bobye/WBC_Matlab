function [ c ] = centroid_sphGD( stride, supp, w, c0)
% Single phase centroid using gradient decent

% Re-prepare
  global A B;
  global stdoutput qpoptim_options;
  
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
  function  [obj grad] = d2energy(x)
  c_supp = reshape(x(1:end-avg_stride),[dim,avg_stride]);
  c_w = x(end-avg_stride+1:end)';
  
  if (any(c_w < 0)) 
      fprintf('%f ',c_w);
      error ('w can not be smaller than zero');  
  end
  grad_supp = zeros(dim, avg_stride);  
  grad_w = zeros(avg_stride,1);
  for it=1:n                              
      [D(it), XX{it}, lambda] = kantorovich(c_supp, c_w, suppx{it}, wx{it});      

      %grad_supp = grad_supp ... 
      %    + repmat(c_w, dim, 1).* c_supp - suppx{it} * XX{it}';

      f=lambda.eqlin(1:avg_stride); f=f-sum(f)/avg_stride;
      grad_w = grad_w + f/sum(abs(f))/(avg_stride^2);
  end
    obj = mean(D);
    grad = [grad_supp(:)*2;grad_w]/n;
    %fprintf(stdoutput, '\n%e', obj );  
  end


  addpath('fminlbfgs_version2c');
  options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs','InitialHessType','identity','GoalsExactAchieve',0);
  [x, ~] = fminlbfgs(@d2energy, [c.supp(:);c.w'], options);
  
  c.supp(:)=x(1:end-avg_stride);
  c.w = x(end-avg_stride+1:end);
end


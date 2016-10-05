function [c] = centroid_sphLP(stride, supp, w, c0, options)
% Single phase centroid using FULL Linear Programming
if isfield(options, 'mosek_path')
    addpath(options.mosek_path);
end
% Re-prepare
  global A B;
  global stdoutput;

  dim = size(supp,1);
  n = length(stride);
  m = length(w);
  if isempty(c0)
    c=centroid_init(stride, supp, w, options);
  else
    c=c0;
  end
  support_size=length(c.w);

  X = zeros(support_size, m);
  D = zeros(n,1);

  iter=0;
  function  obj = d2energy(warm)
  pos=1;
  for it=1:n  
    if warm
    [D(it), X(:, pos:pos+stride(it) -1)] = ... 
    kantorovich(c.supp, c.w, supp(:,pos:pos+stride(it)-1), w(pos:pos+stride(it)-1), ...
    X(:, pos:pos+stride(it) -1));
    else
    [D(it), X(:, pos:pos+stride(it) -1)] = ... 
    kantorovich(c.supp, c.w, supp(:,pos:pos+stride(it)-1), w(pos:pos+stride(it)-1));        
    end
    pos = pos + stride(it);
  end
  obj = mean(D);
  fprintf(stdoutput, '\n\t\t %d\t %e', iter, obj );      
  end

  %d2energy(false);

% optimization

  nIter = 20;
  if isfield(options, 'max_iters') && ~isfield(options, 'support_points')
      nIter=options.max_iters;
  elseif isfield(options, 'support_points')
      nIter = 1;
  end
  suppIter = 1;
  cterm = Inf;
  statusIter = zeros(nIter,1);
  iter_tol = 1E-4;
  for iter=1:nIter
    % update c.supp
    for xsupp=1:suppIter
    d2energy(true);
    if ~isfield(options, 'support_points')
        c.supp = supp * X' ./ repmat(n*c.w, [dim,1]);
    end
    end
    %    x = [ reshape(X, avg_stride*m, 1); c.w'];
    
    % update c.w as well as X, using full LP
    C = pdist2(c.supp', supp', 'sqeuclidean');
    f = reshape(C, support_size*m, 1);
    ff = [f; zeros(support_size,1)];
    
    Aeq = sparse(support_size*n + m + 1, support_size*(m+1));
    beq = sparse(support_size*n + m + 1, 1);

    posi=1;pos=1;posm=1;
    for i=1:n
        stripi = posi:posi+stride(i)-1;
        strips = pos:pos+stride(i)+support_size-1;
        stripsm= posm:posm+stride(i)*support_size-1;
        Aeq(strips,stripsm) = A{support_size, stride(i)};
        beq(strips,1) = [zeros(support_size,1); w(stripi)'];
        Aeq(pos:pos+support_size-1,support_size*m+1:end) = -eye(support_size);
        
        posi= posi + stride(i) ;
        pos = pos + stride(i)+support_size;
        posm = posm + stride(i)*support_size;
    end
    Aeq(pos, posm:end) = ones(1,support_size);
    beq(end) = 1;

    optim_options = optimset('Display','off', 'Diagnostics','off');
    [x, fval, exitflag] = linprog(ff, [], [], Aeq, beq, ...
                                  zeros(support_size*(m+1),1), [], [], ...
                                  optim_options);
    if exitflag < 0
        error('linprog no search direction [%d %f]',exitflag, fval);
    end
    c.w = x(posm:end)'; c.w = c.w/sum(c.w);
    statusIter(iter) = d2energy(false);
    if iter>1 && abs(statusIter(iter)-statusIter(iter-1))<iter_tol*abs(statusIter(iter))
        break;
    end
  end

%  global statusIterRec;
%  statusIterRec(:,2) = statusIter;
  
  %h=figure;
  %plot(statusIter);
  %print(h, '-dpdf', 'centroid_sphLP.pdf');

  %fprintf(stdoutput, ' %f', c.w);
  fprintf(stdoutput, '\n');
  if min(c.w)<1E-8
      warning('There are %d/%d support points whose weights are less than 1E-8\n', sum(c.w<1E-8), support_size);
  end
end
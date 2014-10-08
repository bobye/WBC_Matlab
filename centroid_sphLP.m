function [c] = centroid_sphLP(stride, supp, w)
% Single phase centroid using FULL Linear Programming

% Re-prepare
  global A B;
  global stdoutput;

  dim = size(supp,1);
  n = length(stride);
  m = length(w);

% load start guess
  avg_stride = ceil(mean(stride));
  load(['cstart' num2str(n) '.mat']);

  X = zeros(avg_stride, m);
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

  d2energy(false);

% optimization

  nIter = 20;
  suppIter = 1;
  admmIter = 10;
  cterm = Inf;
  statusIter = zeros(nIter,1);
  for iter=1:nIter
      toc;tic;
    % update c.supp
    for xsupp=1:suppIter
    c.supp = supp * X' ./ repmat(n*c.w, [dim,1]);
    d2energy(true);
    end
    %    x = [ reshape(X, avg_stride*m, 1); c.w'];
    
    % update c.w as well as X, using full LP
    C = pdist2(c.supp', supp', 'sqeuclidean');
    f = reshape(C, avg_stride*m, 1);
    ff = [f; zeros(avg_stride,1)];
    
    Aeq = zeros(avg_stride*n + m + 1, avg_stride*(m+1));
    beq = zeros(avg_stride*n + m + 1, 1);

    posi=1;pos=1;posm=1;
    for i=1:n
        stripi = posi:posi+stride(i)-1;
        strips = pos:pos+stride(i)+avg_stride-1;
        stripsm= posm:posm+stride(i)*avg_stride-1;
        Aeq(strips,stripsm) = A{avg_stride, stride(i)};
        beq(strips,1) = [zeros(avg_stride,1); w(stripi)'];
        Aeq(pos:pos+avg_stride-1,avg_stride*m+1:end) = -eye(avg_stride);
        
        posi= posi + stride(i) ;
        pos = pos + stride(i)+avg_stride;
        posm = posm + stride(i)*avg_stride;
    end
    Aeq(pos, posm:end) = ones(1,avg_stride);
    beq(end) = 1;

    optim_options = optimset('Display','off', 'Diagnostics','off');
    [x, fval, exitflag] = linprog(ff, [], [], Aeq, beq, ...
                                  zeros(avg_stride*(m+1),1), [], [], ...
                                  optim_options);
    if exitflag < 0
        error('linprog no search direction [%d %f]',exitflag, fval);
    end
    c.w = x(posm:end)'; c.w = c.w/sum(c.w);
    statusIter(iter) = d2energy(false);
  end

  global statusIterRec;
  statusIterRec(:,2) = statusIter;
  
  %h=figure;
  %plot(statusIter);
  %print(h, '-dpdf', 'centroid_sphLP.pdf');

  fprintf(stdoutput, ' %f', c.w);
  fprintf(stdoutput, '\n');
end
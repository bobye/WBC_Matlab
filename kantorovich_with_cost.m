function [fval, x, lambda] = kantorovich_with_cost(D, wX, wY, x0)

  global A;
  global default_options optim_options lpoptim_options;
  
  n = length(wX);
  m = length(wY);
  
  wX = wX/sum(wX);
  wY = wY/sum(wY);

  if (length(wX) ~= n || length(wY) ~= m ) 
     error('format not correct');
  end

  % might fail when input is NaN
  if any(isnan(wX)) || any(isnan(wY))
        fprintf('%f ',wX);fprintf('\n');
        fprintf('%f ',wY);fprintf('\n');
  end

%  D = pdist2(X', Y', 'sqeuclidean');
  f = reshape(D, n*m, 1);
  
%  Aeq = A{n,m}(1:end-1,:);beq = [wX'; wY(1:end-1)'];
  Aeq = A{n,m};beq = [wX'; wY'];

  if nargin == 3
      x0 = [];
  else
      x0 = reshape(x0,[n*m,1]);
      if any(isnan(x0))
          fprintf('%f ',x0);fprintf('\n');
          x0=[];
      end
  end

  if any(isnan(f)) || any(isnan(Aeq(:))) || any(isnan(beq))
      disp X;
      disp Y;
      disp D;
      %fprintf('%f ',f);fprintf('\n');
      %fprintf('%f ',Aeq);fprintf('\n');
      %fprintf('%f ',beq);fprintf('\n');
  end
  
  [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, default_options );

  
  if exitflag < 0
  [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, optim_options );
  end
  
  if exitflag < 0
  [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), [], x0, lpoptim_options );
  end
  
  if exitflag < 0
  [x, fval, exitflag, ~, lambda] = linprog(f, [], [], Aeq, beq, zeros(n*m,1), []);
  end

  if exitflag < 0
      save(['err-lingprog-' datestr(clock, 0) '.mat'], 'f', 'Aeq', 'beq', 'n', 'm', 'x0');
      error('linprog no search direction [%d, %f]', exitflag, fval);
  end

  x = reshape(x, n, m);
  x(x<0) = 0;
end

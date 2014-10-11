function [fval, x, lambda] = kantorovich_dual(X, wX, Y, wY)

  global A;
 
  n = size(X,2);
  m = size(Y,2);

  wX = wX/sum(wX);
  wY = wY/sum(wY);

  if (length(wX) ~= n || length(wY) ~= m ) 
     error('format not correct');
  end

  D = pdist2(X', Y', 'sqeuclidean');
  b = reshape(D, n*m, 1);
  
  Aineq = A{n,m}'; f = -[wX'; wY'];
  Aeq = [ones(1,n),zeros(1,m)]; beq = 0;
  
  [x, fval, ~, ~, lambda] = linprog(f, Aineq, b, Aeq, beq );

end

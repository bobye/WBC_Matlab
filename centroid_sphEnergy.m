function [ c ] = centroid_sphEnergy( stride, supp, w, c0 )
% Compute total cost only and return c0
  n = length(stride);
  posvec=[1,cumsum(stride)+1];
  c=c0;
  D = zeros(n,1);  
  for it=1:n                              
    idx=posvec(it):posvec(it+1)-1;  
    [D(it), ~] = kantorovich(c.supp, c.w, supp(:,idx), w(idx));
  end
  obj = mean(D);
  fprintf('\n\t\t %f', obj ); 

end


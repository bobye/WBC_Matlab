function [model, beta] = mixture_sphSimulAnn(stride, supp, w, ncomp)

  d = size(supp,1);
  n = length(stride); % number of instances
  m = length(w);
  posvec=[1,cumsum(stride)+1];
  avg_stride = ceil(mean(stride));  
 
  %% Wasserstein barycenter initial
  if (ncomp == 1)
  load(['cstart' num2str(n) '.mat']);  
  model.ncomp = ncomp;
  model.supp = c.supp;
  model.w = c.w;  
  beta=ones(ncomp, n)/ncomp;
  %% Wasserstein mixture initial
  else
  load(['cstart' num2str(n) '.mat']);
  model.ncomp = ncomp;
  model.supp = c.supp;
  for i=1:ncomp
    model.w=gamrnd(100*ones(ncomp, avg_stride), 1.);
    model.w=bsxfun(@times, model.w, 1./sum(model.w, 2));
    beta=ones(ncomp, n)/ncomp;
  end
  end
  %%
  assert(ncomp == size(model.w,1));
  

  nIter = 3000;
  sigma=power(1/nIter, 1/nIter);
  
  suppx{n,1}=[];
  wx{n,1}=[];
  strips{n,1}=[];
  f1{n,1}=[];
  f2{n,1}=[];
  U{n,1}=[];
  L{n,1}=[];
  for iter=1:n
      strips{iter} = posvec(iter):(posvec(iter)+stride(iter)-1);
      suppx{iter} = supp(:,strips{iter});
      wx{iter} = w(strips{iter});
      f1{iter} = zeros(avg_stride,1);
  end
  
  C = pdist2(model.supp', supp', 'sqeuclidean'); mC=mean(C(:));
  gamma = 2E-3;
  T=1E-1 * mC;
  obj=zeros(nIter,1);
  pobj=zeros(nIter,1);
  maxobj=-Inf;
  md=zeros(avg_stride, 1);
  msupp=zeros(d, avg_stride);
  for iter=1:nIter
      order=randperm(n);
      csupp=zeros(d, avg_stride);
      csuppw=zeros(1, avg_stride);
      for ii=1:n
          i=order(ii); % shuffling samples every epoch
          % Gibbs sampler
          W=model.w' * beta(:,i);
          Ci=C(:,strips{i});
          Cl=bsxfun(@plus, -Ci, f1{i});
          [L{i}, midx2] = max(Cl, [], 1);
          f2{i} = L{i} + exprnd(T./wx{i} , [1, stride(i)]);
          Cu=bsxfun(@plus,  Ci, f2{i});
          [U{i}, midx1] = min(Cu, [], 2);
          f1{i} = U{i} - exprnd(T./W, [avg_stride, 1]);
          % mirror descent          
          obj(iter)=obj(iter) + dot(f1{i}, W) - dot(f2{i}, wx{i});          
          gd=f1{i} - mean(f2{i}); % subgradient
          md=0.9*md + gd;
          if maxobj > 0          
          model.w = model.w .* exp(-(gamma/mC) * beta(:,i) * md');
          model.w = bsxfun(@times, model.w, 1./(sum(model.w, 2)+eps));
          beta(:,i) = beta(:,i) .* exp(-(gamma/mC) * (model.w * md) );
          beta(:,i) = beta(:,i) / (sum(beta(:,i)) + eps); 
          end
                   
          
          % support points
          % method 1
          if true
          Pi= zeros(avg_stride,stride(i));
          midx1=sub2ind([avg_stride,stride(i)], 1:avg_stride, midx1');
          midx2=sub2ind([avg_stride,stride(i)], midx2, 1:stride(i));
          Pi(midx1)=W; Pi(midx2)=Pi(midx2)+wx{i}; Pi=Pi/2;
          end
          % method 2
          if false
          s2=exp(-Cu/T); s2=bsxfun(@times, s2, 1./sum(s2,2));
          s1=exp(-Cl/T); s1=bsxfun(@times, s1, 1./sum(s1,1));
          s2=bsxfun(@times, s2, W);
          s1=bsxfun(@times, s1, wx{i});
          Pi=(s1+s2) / (sum(s1(:)) + sum(s2(:)));
          end
          pobj(iter) = pobj(iter) + dot(Pi(:), Ci(:));
          csupp = csupp + suppx{i} * Pi';          
          csuppw = csuppw + sum(Pi, 2)';
      end
      obj(iter)=obj(iter)/n;
      pobj(iter)=pobj(iter)/n;
      fprintf('OBJ: %f ', obj(iter));
      fprintf('POBJ: %f ', pobj(iter));
      %fprintf('W: '); fprintf('%f ', model.w); 
      if (mod(iter, 10) == 0 || (obj(iter) > maxobj)) && obj(iter) > 0
      csupp = bsxfun(@times, csupp, 1./csuppw);
      msupp = 0.9 * msupp + (csupp - model.supp);
      model.supp = model.supp + gamma * msupp;
      C = pdist2(model.supp', supp', 'sqeuclidean'); mC=mean(C(:));
      maxobj=obj(iter);
      fprintf('SU');
      end
      %model.supp     
      fprintf('\n');      
      if (obj(iter) < pobj(iter))
        T=T*sigma;      
      else
        T=T*sigma*sigma;
        gamma=gamma*sigma*sigma;
      end
  end
  plot(obj); hold on; plot(pobj);
  [obj(end),pobj(end)]
end
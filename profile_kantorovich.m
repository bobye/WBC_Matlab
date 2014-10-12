%% Profile the speeds of W distance
%% Load data
clear;
setparam;
size=20;

s_modalities = 2;
d_modalities = [3, 3];
filename='../total.txt';
db = loaddata(size, s_modalities, d_modalities, filename);

max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(max_stride);

% %% 
tic;
  for s=1:s_modalities
      idx=[1,cumsum(db{s}.stride)];
      for i=1:size
          for j=i+1:size
              d1.supp = db{s}.supp(:,idx(i):idx(i+1)-1);
              d1.w = db{s}.w(idx(i):idx(i+1)-1);
              d2.supp = db{s}.supp(:,idx(j):idx(j+1)-1);
              d2.w = db{s}.w(idx(j):idx(j+1)-1);
              [f0,x0,lambda0] = kantorovich(d1.supp, d1.w, d2.supp, d2.w);
              %[f1,x1,lambda1] = kantorovich_dual(d1.supp, d1.w, d2.supp, d2.w);
          end
      end
  end
toc
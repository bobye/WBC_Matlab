%% Precomputation
function kantorovich_prepare(max_stride1, max_stride2)

global A B;

A = cell(max_stride1);
B = cell(max_stride2);

for n=1:max_stride1
  for m=1:max_stride2
    A{n,m} = sparse(n+m, n*m);
    for k=1:n
      A{n,m}(k, k:n:n*m) = 1;
    end
    for k=1:m
      A{n,m}(n+k, (k-1)*n +(1:n)) = 1;
    end
  end
end


for n=1:max_stride1
  for m=1:max_stride2
    B{n,m} = sparse(kron(ones(m), eye(n)));
  end
end

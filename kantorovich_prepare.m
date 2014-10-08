%% Precomputation
function kantorovich_prepare(max_stride)

global A B;

for n=1:max_stride
  for m=1:max_stride
    A{n,m} = zeros(n+m, n*m);
    for k=1:n
      A{n,m}(k, k:n:n*m) = 1;
    end
    for k=1:m
      A{n,m}(n+k, (k-1)*n +(1:n)) = 1;
    end
  end
end


for n=1:max_stride
  for m=1:max_stride
    B{n,m} = kron(ones(m), eye(n));
  end
end

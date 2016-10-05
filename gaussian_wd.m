function D = gaussian_wd(V1, V2)
% size(V1) = [d*d, n1]
% size(V2) = [d*d, n2]
% size(D) = [n1, n2]


d=sqrt(size(V1,1));
assert(size(V1,1) == d*d);
assert(size(V2,1) == d*d);
n1=size(V1,2);
n2=size(V2,2);

D=zeros(n1, n2);
V1=reshape(V1, [d d n1]);
V2=reshape(V2, [d d n2]);
for i=1:n1
    d1=sqrtm(V1(:,:,i));
    t1=trace(V1(:,:,i));
    for j=1:n2
        v2=V2(:,:,j);
        d2=sqrtm(d1 * v2 * d1);
        D(i,j) = t1 + trace(v2) - 2 * trace(d2);
    end
end
end
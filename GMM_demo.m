N=100;
m=3;
d=2;
Tau = [1 .5; .5 2];
df =10;

mu = [2,3];
sigma = [1,1.5;1.5,3];

supp=zeros(d+d*d, m*N);
w = zeros(1,m*N);
for i=1:N
    for j=1:m
        r = mvnrnd(mu,sigma,1)';
        S = iwishrnd(Tau,df)*(df-2-1);
        supp(:,(i-1)*m+j)=[r;S(:)];       
    end
    t=rand(m,1);
    w((i-1)*m+1:i*m) = t/sum(t);    
end
stride=m*ones(1,N);

%% 
c0.supp=supp(:, 1:m);
c0.w=w(1:m);
options.badmm_max_iters=2000;
c=centroid_sphBregman_GMM(stride, supp, w, c0, options);


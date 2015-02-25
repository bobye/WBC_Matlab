function [] = syntheticdata(dim, stride, meta_class, total_size)
% generate synthetic data

centroids.m = stride;
centroids.num = meta_class;
samples.num = total_size;


perturbation.supp = .1/(centroids.m^(1/dim));
perturbation.w = .5/centroids.num;

centroids.supp = randn(dim, centroids.m, centroids.num);
centroids.w = gamrnd(ones(centroids.m,centroids.num), 1); 
centroids.w = centroids.w ./ repmat(sum(centroids.w), [centroids.m, 1]);


samples.labels = randi(centroids.num, [samples.num, 1]);

samples.supp = reshape(centroids.supp(:,:,samples.labels), dim, centroids.m*samples.num);
samples.w = centroids.w(:,samples.labels);

samples.supp = samples.supp + perturbation.supp * trnd(2*ones(dim, centroids.m*samples.num));
samples.w = samples.w + perturbation.w * randn([centroids.m, samples.num]);
samples.w(samples.w<(1E-3 / centroids.m)) = 1E-3 / centroids.m;
samples.w = samples.w ./ repmat(sum(samples.w), [centroids.m, 1]);

filename = ['../data/synthetic_data/' num2str(samples.num) '_' ... 
             num2str(dim) '_' ...
             num2str(centroids.m) '_' ...
             num2str(centroids.num)];
%fid = fopen([filename '.d2'],'w');
for i=1:samples.num
    printf( '%d\n', dim);
    printf( '%d\n', centroids.m);
    printf( '%f ', samples.w(:,i)'); printf('\n');
    for j=(centroids.m*(i-1)+1) : centroids.m*i
        printf( '%f ',samples.supp(:, j)); printf( '\n');
    end
    printf('\n');
end
         
%fclose(fid);

%fid = fopen([filename '_c.d2'],'w');
%for i=1:centroids.num
    %printf( '%d\n', dim);
    %printf( '%d\n', centroids.m);
    %printf( '%f ', centroids.w(:,i)'); printf( '\n');
    %for j=1:centroids.m
        %printf( '%f ', centroids.supp(:,j,i)); printf( '\n');
    %end
%end
%fclose(fid);



% generate synthetic data

dim = 3;
centroids.m = 5;
centroids.num = 10;
samples.num = 3000;


perturbation.supp = .1/centroids.m^(1/dim);
perturbation.w = .5/centroids.num;

centroids.supp = rand(dim, centroids.m, centroids.num);
centroids.w = rand(centroids.m, centroids.num); 
centroids.w = centroids.w ./ repmat(sum(centroids.w), [centroids.m, 1]);


samples.labels = randi(centroids.num, [samples.num, 1]);

samples.supp = reshape(centroids.supp(:,:,samples.labels), dim, centroids.m*samples.num);
samples.w = centroids.w(:,samples.labels);

samples.supp = samples.supp + perturbation.supp * randn([dim, centroids.m*samples.num]);
samples.w = samples.w + perturbation.w * randn([centroids.m, samples.num]);
samples.w(samples.w<0) = 0;
samples.w = samples.w ./ repmat(sum(samples.w), [centroids.m, 1]);

filename = ['../' num2str(samples.num) '_' ... 
             num2str(dim) '_' ...
             num2str(centroids.m) '_' ...
             num2str(centroids.num)];
fid = fopen([filename '.txt'],'w');
for i=1:samples.num
    fprintf(fid, '%d\n', dim);
    fprintf(fid, '%d\n', centroids.m);
    fprintf(fid, '%f ', samples.w(:,i)'); fprintf(fid,'\n');
    for j=(centroids.m*(i-1)+1) : centroids.m*i
        fprintf(fid, '%f ',samples.supp(:, j)); fprintf(fid, '\n');
    end
    fprintf(fid,'\n');
end
         
fclose(fid);

fid = fopen([filename '_c.txt'],'w');
for i=1:centroids.num
    fprintf(fid, '%d\n', dim);
    fprintf(fid, '%d\n', centroids.m);
    fprintf(fid, '%f ', centroids.w(:,i)'); fprintf(fid, '\n');
    for j=1:centroids.m
        fprintf(fid, '%f ', centroids.supp(:,j,i)); fprintf(fid, '\n');
    end
end
fclose(fid);



%clear;
numOfSamples = 5000;
if ~exist('numOfSamples')
    numOfSamples = 50;
end
setparam;
%% Load data

fprintf(stdoutput, 'Loading data ... ');

s_modalities = 2;
d_modalities = [3, 3];
filename='../total.txt';
db = loaddata(numOfSamples, s_modalities, d_modalities, filename);

%%
global statusIterRec;

max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(max_stride);

matlabpool('open', num_of_cores); % start parallel workers
clusters = d2clusters(db, 6);
matlabpool('close');

save clusters.dat clusters IDX labels

n = size(statusIterRec,1);

%h = figure;
%plot((1:n)', statusIterRec(:,1),'-or', ...
%     (1:n)', statusIterRec(:,2),'-dg', ...
%     (1:n)', statusIterRec(:,3),'-+b');
 
%plot((1:n)', statusIterRec(:,1),'-or', ...
%     (1:n)', statusIterRec(:,2),'-dg');
 
%err = sqrt(kantorovich(bufferc{1}.supp, bufferc{1}.w, bufferc{2}.supp, bufferc{2}.w)) ... 
%    /norm(bufferc{2}.supp,'fro');

%print(h, '-dpdf', ['centroid_sphALL' num2str(numOfSamples) '.pdf']);

fprintf('%d %d %f %f', numOfSamples, num_of_cores, ctime(1), ctime(2));
%err


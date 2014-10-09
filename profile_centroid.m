%% Profile the speed of W distance
%% Load data
clear;
setparam;
size=100;

s_modalities = 2;
d_modalities = [3, 3];
filename='../total.txt';
db = loaddata(size, s_modalities, d_modalities, filename);

max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(max_stride);

% %% 
% tic;
% for s=1:s_modalities
%     idx=[1,cumsum(db{s}.stride)];
%     for i=1:size
%         for j=i+1:size
%             [~,x0] = kantorovich(db{s}.supp(:,idx(i):idx(i+1)), ...
%                                  db{s}.w(idx(i):idx(i+1)), ...
%                                  db{s}.supp(:,idx(j):idx(j+1)), ...
%                                  db{s}.w(idx(j):idx(j+1)));                         
%         end
%     end
% end
% toc

%%
for s=1:s_modalities
    centroid_singlephase(db{s}.stride, db{s}.supp, db{s}.w, []);
end

%%
for s=1:s_modalities
    centroid_sphLP(db{s}.stride, db{s}.supp, db{s}.w);
end

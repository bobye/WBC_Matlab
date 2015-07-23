%% Profile the speeds of W distance
%% Load data
clear;
setparam;
size=1000;

s_modalities = 2;
d_modalities = [3 3];
filename='../data/mountaindat.d2';
db = loaddata(size, s_modalities, d_modalities, filename);

max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(max_stride);


%%
%profile on
%for s=1:1 %1:s_modalities
%    c.GD = centroid_sphGD(db{s}.stride, db{s}.supp, db{s}.w, []);
%end
%profile off
%profile viewer
%%
profile on
for s=1:1 %1:s_modalities
    c.Bregman = centroid_sphBregman(db{s}.stride, db{s}.supp, db{s}.w, []);
    centroid_sphEnergy(db{s}.stride, db{s}.supp, db{s}.w, c.Bregman);    
end
profile off
profile viewer

%%
%profile on
%for s=2:2 %1:s_modalities
%    c.ADMM = centroid_sphADMM(db{s}.stride, db{s}.supp, db{s}.w, []);
%end
%profile off
%profile viewer
%%
%for s=1:1 %1:s_modalities
%    c.LP = centroid_sphLP(db{s}.stride, db{s}.supp, db{s}.w);
%end

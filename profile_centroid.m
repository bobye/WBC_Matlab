%% Profile the speeds of W distance
%% Load data
clear;
setparam;
max_sample_size=1000;

s_modalities = 1; %2;
d_modalities = 2; %[3 3];
filename='weather45.d2';
db = loaddata(max_sample_size, s_modalities, d_modalities, filename);

%% Additional Configuration
options.max_support_size=50;
options.mosek_path='/Users/jxy198/mosek/7/toolbox/r2013a/'; % optional
options.init_method='kmeans';
options.support_size=50;
options.method='badmm'; % {'lp', 'gd', 'badmm', 'admm', 'ibp'}

%% Compute Wasserstein Barycenter
[c, OT]=Wasserstein_Barycenter(db, [], options);

%% write results
write_centroid_and_ot(c, OT, 'weather45');


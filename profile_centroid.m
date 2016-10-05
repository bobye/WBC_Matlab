%% Profile the speeds of W distance
%% Load data
clear;
setparam;
max_sample_size=1000;

s_modalities = 2;
d_modalities = [3 3];
filename='mountaindat.d2';
db_tmp = loaddata(max_sample_size, s_modalities, d_modalities, filename);


db=cell(1,1); db{1}=db_tmp{1}; % only use the color quantization data

%% Additional Configuration
options.max_support_size=50;
options.mosek_path='/Users/jxy198/mosek/7/toolbox/r2013a/'; % optional
options.init_method='kmeans'; % {'kmeans', 'mvnrnd'}
options.support_size=60;
options.max_support_size=options.support_size;

%% Set Initialization
n=length(db);
c0=cell(n,1);
for s=1:n
    c0{s}=centroid_init(db{s}.stride, db{s}.supp, db{s}.w, options);
end


%% comment out this line to use fixed support points of barycenter
%% options.support_points= 'fixed';

%% Compute Wasserstein Barycenter (LP)

%options.method='lp'; % {'lp', 'gd', 'badmm', 'admm', 'ibp'}
%[c, OT]=Wasserstein_Barycenter(db, c0, options);


%% Compute Wasserstein Barycenter (B-ADMM)
options.method='badmm'; % {'lp', 'gd', 'badmm', 'admm', 'ibp'}
options.badmm_tol=1e-3;
[c, OT]=Wasserstein_Barycenter(db, c0, options);

%% Compute Wasserstein Barycenter (IBP)
options.method='ibp'; % {'lp', 'gd', 'badmm', 'admm', 'ibp'}
for vareps = [0.2, 0.1, 0.02, 0.01, 0.005]
    options.ibp_vareps=vareps;
    options.ibp_max_iters=1e6;
    options.ibp_tol=1e-3; % set ibp_tol higher to have faster performance
    [c, OT]=Wasserstein_Barycenter(db, c0, options);
end
%% write results
%write_centroid_and_ot(c, OT, 'weather45');


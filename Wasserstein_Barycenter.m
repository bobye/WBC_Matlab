function [c, OT] = Wasserstein_Barycenter(db, c0, options)
% Compute Wasserstein barycenter with different methods
% Input:
%  db -- a struct that stores n-phase discrete distribution data.
%           db{1...n} are data associated to different phases;
%           db{i}.stride is an array of the support sizes of individual
%           instance at phase i;
%           db{i}.supp is a matrix storing all support points across
%           instances by colums at phase i;
%           db{i}.w is an array storing all weights across instances at
%           phase i;
%           It is required that length(db{i}.w)==size(db{i}.supp,2).
%  c0 -- a multi-phase discrete distribution as the initial start
%           c0.supp is a matrix storing support points by columns
%           c0.w is an array of weights. It can be automatically generated
%           if it is set to []. 
%           
%  options -- 
%          shared:
%           options.support_size: the desired support size of barycenter
%
%           options.support_points: if specified, the support of barycenter 
%           is fixed to options.support_points. Then it is required that
%           options.support_size == size(options.support_points, 2)
% 
%           options.init_method: {'kmeans', 'mvnrnd'} There are two options
%           for initialization of barycenter, if c0 is empty. One is
%           Kmeans, and the other is multivariate normal.
%
%           options.max_support_size: the maximum desired support size of
%           barycenter, typically set to 3*support_size(c).
%
%           options.method: {'lp', 'gd', 'badmm' (default), 'admm', 'ibp'}
%
%          method specific options:
%           'lp':
%             options.max_iters (default, 20)
%           'gd':
%             options.max_iters (default, 20)
%             options.gd_step_size (default, 0.008)
%           'badmm':
%             options.badmm_max_iters (default, 2000)
%             options.badmm_rho (default, 2.0)
%             options.badmm_tau (default, 10)
%           'admm':
%             options.admm_max_iters (default, 50)
%             options.admm_inner_iters (default, 10)
%             options.admm_rho (default, 50.)
%           'ibp':
%             options.ibp_max_iters (default, 20000)
%             options.ibp_vareps (default, 0.01)
%             options.ibp_tol (default, 1E-4)
%
% Output:
%  c  -- Wasserstein barycenter
%  OT -- matching matrix between c and each instance




max_stride = max(cellfun(@(x) max(x.stride), db));
kantorovich_prepare(options.max_support_size,max_stride);

method='badmm';
if isfield(options, 'method')
    method = options.method;
end

n=length(db);
c=cell(n,1);
OT=cell(n,1);
for s=1:n
    if strcmp(method, 'lp')
        c{s}=centroid_sphLP(db{s}.stride, db{s}.supp, db{s}.w, options);
    elseif strcmp(method, 'gd')
        c{s}=centroid_sphGD(db{s}.stride, db{s}.supp, db{s}.w, c0, options);
    elseif strcmp(method, 'badmm')
        c{s}=centroid_sphBregman(db{s}.stride, db{s}.supp, db{s}.w, c0, options);
    elseif strcmp(method, 'admm')
        c{s}=centroid_sphADMM(db{s}.stride, db{s}.supp, db{s}.w, c0, options);
    elseif strcmp(method, 'ibp')
        c{s}=centroid_sphIBP(db{s}.stride, db{s}.supp, db{s}.w, c0, options);
    end
    [~, OT{s}] =centroid_sphEnergy(db{s}.stride, db{s}.supp, db{s}.w, c{s});
end
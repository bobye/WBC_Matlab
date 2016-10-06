%% Set global param
global stdoutput IDX ctime default_options optim_options lpoptim_options qpoptim_options bufferc num_of_cores; 
stdoutput = 1;
ctime=zeros(2,1);
num_of_cores = 12;
optim_options   = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off');
lpoptim_options = optimoptions('linprog','Algorithm','dual-simplex','Display','off', 'Diagnostics','off');
%qpoptim_options = optimset('Display','off', 'LargeScale','off', 'Diagnostics','off', 'Algorithm','active-set');
qpoptim_options = optimset('Display','off', 'Diagnostics','off', 'Algorithm','interior-point-convex');


default_options = optimset('Display','off', 'Diagnostics','off');

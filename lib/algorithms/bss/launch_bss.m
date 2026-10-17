% LAUNCH_BSS runs the BSS baseline (contrib-bss) on a QSI test case
%
% SYNTAX: launch_bss (PROBLEM, CONFIG, IT, DATA_DIR)
%
%   runs Bayesian Subset Simulation on the test case PROBLEM, starting from
%   the initial design of run IT in DATA_DIR/doe_init, and writes the final
%   design to DATA_DIR/results/design/doe_BSS_<name>_<IT>_1.csv.

function launch_bss (problem, config, it, data_dir)

problem_bss = convert_problem_to_bss_format (problem);

% Construct parameter estimation method
[prm, f, ~] = problem ();
param_estim = @(not_used, xi, yi) ...
    estim_matern (xi, yi, prm.list_cov, config.lognugget);

% BSS options
options = bss_get_default_options ([], problem_bss);
options.param_estim = param_estim;
options.samplesize = config.bss_options.samplesize;
options.verbose = config.bss_options.verbose;

fprintf ('Run %d\n', it);

doe_init = fullfile (data_dir, "doe_init", ...
    sprintf ("doe_init_%s_%d_init.csv", prm.name, it));
x0 = stk_dataframe (readmatrix_ (doe_init));
y0 = f (x0);

[~, doe] = bss (problem_bss, options, x0, y0);

doe_file = fullfile (data_dir, "results/design", ...
    sprintf("doe_BSS_%s_%d_1.csv", prm.name, it));
writematrix_ (double (doe.xi), doe_file)

end % function


function problem_bss = convert_problem_to_bss_format (problem)

[prm, f, ~] = problem ();

% Input distribution: uniform on the hyper-rectangle
% NB: We use a uniform distribution on U = X x S, regardless
%     of the actual distribution of "s" in the problem
a = double (prm.BOX(1, :));
b = double (prm.BOX(2, :));
input_distrib = UniformDistribution (a, b);

% Target threshold
u_final = prm.const(2);
assert (isequal (prm.const, [-inf; u_final]));

problem_bss.name = prm.name;
problem_bss.dim = prm.dim_x + prm.dim_s;
problem_bss.cost_fun = f;
problem_bss.input_distrib = input_distrib;
problem_bss.u_final = u_final;
problem_bss.Pf_true = 0;

end % function

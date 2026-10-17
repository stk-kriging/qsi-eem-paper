% Cost of the QSI-SUR and QSI-EEM criteria against the grid size
% M = |Xt x St|, with no importance sampling and no candidate ranking.
% Writes ./outputs/results/time/criterion_time.csv; the figure is drawn by
% make_figure_time_criterion_eem_sur.
%
% QSI-SUR peaks near 14 GB at M = 25000, chol and the regularization loop
% of stk_cholcov each duplicating the covariance matrix; a grid size that
% fails is recorded as NaN and the sweep continues.  Watch for the
% stk_cholcov warning if the test case or the nugget changes: those
% retries would be counted in the QSI-SUR timings.

close all;

fprintf ('Launching QSI-SUR / QSI-EEM criterion runtime comparison\n');

script_dir = fileparts (mfilename ('fullpath'));
paper_code_dir = fileparts (script_dir);

data_dir = fullfile (script_dir, 'outputs');
time_dir = fullfile (data_dir, 'results', 'time');

if ~ isfolder (time_dir)
    mkdir (time_dir);
end

debug_mode = false;

case_name = 'branin_mod';   % f_1 of Section 4.3
run_id = 1;                 % initial design used to fit the GP

% pts_s points in S times a varying number in X: grid size 1000 ... 25000
pts_s = 100;
pts_x_list = 10:10:250;
nb_repetitions = 10;

if debug_mode
    pts_x_list = 10:10:30;
    nb_repetitions = 2;
end

grid_sizes = pts_x_list * pts_s;

[problem, config] = paper_get_test_case (case_name, false);
[prm, f, s_trnsf] = problem ();

config.pts_s = pts_s;

fprintf ('%s: grid sizes %d ... %d (%d values), %d repetition(s)\n', ...
    case_name, grid_sizes(1), grid_sizes(end), numel (grid_sizes), ...
    nb_repetitions);

% Gaussian quadrature used by the QSI-SUR criterion (scalar output)
if prm.M ~= 1
    error ('This benchmark assumes a scalar output (prm.M == 1).');
end
quantOpt.nbLevels = config.nVar;
quantOpt.useGaussHermite = 1;

% design and model fitted once: this compares criteria, not estimation
paper_set_seed (run_id);
paper_generate_initial_designs (problem, config, run_id, data_dir);

file_grid = sprintf ('doe_init_%s_%d_init.csv', prm.name, run_id);
di = readmatrix_ (fullfile (data_dir, 'doe_init', file_grid));
zi = f (di);

dn = stk_dataframe (di);
zn = stk_dataframe (zi);

Model = stk_model ();
for m = 1:prm.M
    Model(m) = estim_matern (dn, zn(:,m), prm.list_cov, config.lognugget);
end

fprintf ('Initial design: %d points, %d + %d dimensions\n', ...
    size (dn, 1), prm.dim_x, prm.dim_s);

% warm-up: the first call pays JIT and the quadrature cache of modif_mu
warmup_x = pts_x_list(1);
config.keep_x = warmup_x;
xt = stk_sampling_randomlhs (warmup_x, prm.dim_x, prm.BOXx);
st = s_trnsf (stk_sampling_randomlhs (pts_s, prm.dim_s, prm.BOXs));
dt = adapt_set (xt, st);
cand_index = max (1, round (size (dt, 1) / 2));
time_sur (Model, dn, zn, dt, dt(cand_index,:), cand_index, ...
    ones (warmup_x, config.nVar), prm, config, quantOpt);
time_eem (Model, dn, zn, dt, dt(cand_index,:), ...
    ones (warmup_x, 1) / warmup_x, prm, config);

results = nan (numel (grid_sizes) * nb_repetitions, 6);
row = 0;

t_all = tic;

for g = 1:numel (grid_sizes)

    pts_x = pts_x_list(g);
    M = grid_sizes(g);

    % keep_x is read by the reshape operations of the QSI-SUR criterion
    config.keep_x = pts_x;

    fprintf ('\n=== grid size %d (%d x %d), %d/%d\n', ...
        M, pts_x, pts_s, g, numel (grid_sizes));

    for rep = 1:nb_repetitions

        row = row + 1;
        paper_set_seed (1000 * g + rep);

        % Plain grid in X x S: no importance sampling, no ranking
        xt = stk_sampling_randomlhs (pts_x, prm.dim_x, prm.BOXx);
        st = stk_sampling_randomlhs (pts_s, prm.dim_s, prm.BOXs);
        st = s_trnsf (st);
        dt = adapt_set (xt, st);

        % Candidate point: a fixed point of the grid, away from its ends
        cand_index = max (1, round (size (dt, 1) / 2));
        pt = dt(cand_index, :);

        % Uniform integration weights
        IS_sur = ones (pts_x, config.nVar);
        IS_eem = ones (pts_x, 1) / pts_x;

        [t_sur_prep, t_sur_eval] = time_sur ...
            (Model, dn, zn, dt, pt, cand_index, IS_sur, prm, config, quantOpt);

        [t_eem_prep, t_eem_eval] = time_eem ...
            (Model, dn, zn, dt, pt, IS_eem, prm, config);

        results(row, :) = [M, rep, t_sur_prep, t_sur_eval, ...
            t_eem_prep, t_eem_eval];

        fprintf ('  rep %d | QSI-SUR %8.3f s (prep %.3f) | QSI-EEM %8.4f s (prep %.4f)\n', ...
            rep, t_sur_prep + t_sur_eval, t_sur_prep, ...
            t_eem_prep + t_eem_eval, t_eem_prep);

    end

end

fprintf ('\n=== done in %.1f min\n', toc (t_all) / 60);

out_file = fullfile (time_dir, 'criterion_time.csv');
writematrix_ (results, out_file);
fprintf ('Wrote %s\n', out_file);

manifest = cell (5, 1);
manifest{1} = sprintf ('debug_mode=%d', debug_mode);
manifest{2} = sprintf ('case=%s', case_name);
manifest{3} = sprintf ('grid_sizes=%s', mat2str (grid_sizes));
manifest{4} = sprintf ('nb_repetitions=%d', nb_repetitions);
manifest{5} = sprintf ('nTraj=%d, nVar=%d', config.nTraj, config.nVar);
paper_write_manifest (data_dir, manifest);

cd (paper_code_dir);


function [t_prepare, t_evaluate] = time_sur ...
    (Model, dn, zn, dt, pt, cand_index, IS, prm, config, quantOpt)

t_prepare = NaN;
t_evaluate = NaN;

n = size (dn, 1);
nb_dt = size (dt, 1);

try

    % preparation: sample paths over [dn; dt], i.e. the Cholesky factorization
    t0 = tic;
    traj = zeros (n + nb_dt, config.nTraj, prm.M);
    for m = 1:prm.M
        traj(:,:,m) = stk_generate_samplepaths ...
            (Model(m), [dn; dt], config.nTraj);
    end
    t_prepare = toc (t0);

    % evaluation at one candidate point
    t0 = tic;

    var = [];
    for m = 1:prm.M
        p = stk_predict (Model(m), dn, zn(:,m), pt);
        [var, weight] = quantization (p.mean, sqrt (p.var), quantOpt);
    end

    xc_ind = n;
    bool = cell (config.nVar, prm.M);

    for m = 1:prm.M

        [~, lambda] = stk_predict (Model(m), [dn; pt], [], [dn; dt]);
        lambda_dn = lambda(1:(size (lambda, 1) - 1), :);
        lambda_pt = lambda(size (lambda, 1), (xc_ind + 1):size (traj, 1));

        trajCond_dn = stk_conditioning ...
            (lambda_dn, zn(:,m), traj(:,:,m), 1:xc_ind);
        trajCond_dn = trajCond_dn((xc_ind + 1):size (trajCond_dn, 1), :);

        lambda_pt = reshape (lambda_pt, config.keep_x, config.pts_s);
        tensor_dn = reshape ...
            (trajCond_dn, config.keep_x, config.pts_s, config.nTraj);

        traj_dt = traj((xc_ind + 1):size (traj, 1), :, m);
        traj_ind = reshape (traj_dt(cand_index, :), 1, 1, config.nTraj);

        for k = 1:config.nVar
            bool{k,m} = check_constraints_trajs ...
                (tensor_dn, traj_ind, var(k), lambda_pt, prm.const(:,m));
        end

    end

    proba = proba_tau_trajs (bool, prm.alpha);
    s = min (proba, 1 - proba);          % config.critName == "m"
    crit = sum (weight .* (mean (IS .* s, 1)));  %#ok<NASGU>

    t_evaluate = toc (t0);

catch err
    warning ('QSI-SUR at grid size %d failed: %s', nb_dt, err.message);
end

end % function


function [t_prepare, t_evaluate] = time_eem ...
    (Model, dn, zn, dt, pt, IS, prm, config)

t_prepare = NaN;
t_evaluate = NaN;

nb_dt = size (dt, 1);

try

    % preparation: posterior model, kriging weights and plug-in estimator
    t0 = tic;

    precomp = struct ();
    for m = 1:prm.M
        precomp.post_model{m} = stk_model_gpposterior (Model(m), dn, zn(:,m));
        kreq_dt = stk_make_kreq (precomp.post_model{m}, dt);
        precomp.lambda_mu_dt{m} = (get (kreq_dt, 'lambda_mu'))';
    end

    mu_n = predict_mean_var (Model, dn, zn, dt, prm);
    estim_current = get_estim (mu_n, prm.const, prm, config);

    t_prepare = toc (t0);

    % evaluation at one candidate point
    t0 = tic;

    crit = modif_mu (pt, [], Model, dn, zn, dt, estim_current, mu_n, ...
        IS, prm.const, prm, config, precomp);  %#ok<NASGU>

    t_evaluate = toc (t0);

catch err
    warning ('QSI-EEM at grid size %d failed: %s', nb_dt, err.message);
end

end % function

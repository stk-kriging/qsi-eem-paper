% Criterion maps for the QSI-SUR / QSI-EEM comparison.  Along one QSI-EEM
% run on branin_mod, both criteria are evaluated at the same posteriors on
% the same regular candidate grid.  Writes the maps as EPS files, and a CSV
% with the Spearman correlation between the maps and the location of both
% optima.

close all;

script_dir = fileparts (mfilename ('fullpath'));

data_dir    = fullfile (script_dir, 'outputs');
results_dir = fullfile (data_dir, 'results');
fig_dir     = fullfile (results_dir, 'figures');
stats_dir   = fullfile (results_dir, 'stats');
paper_prepare_data_dir (data_dir);
if ~ isfolder (stats_dir)
    mkdir (stats_dir);
end

case_name = 'branin_mod';
it = 10;                    % run id, also the seed
grid_pts = 50;              % candidate grid of grid_pts^2 points
compare_steps = [0 5 10];   % steps (after the initial design) with maps
PTS_REF = 250;              % reference grid for the overlays

opts = struct ( ...
    'fig_width_cm',  8.0, ...   % included at 70mm in the paper
    'fig_height_cm', 6.4, ...
    'font_size',     9);

[problem, config] = paper_get_test_case (case_name);
[prm, f, s_trnsf] = problem ();
config.batch_size = 1;

paper_generate_initial_designs (problem, config, it, data_dir);
paper_set_seed (it);

file_grid = sprintf ('doe_init_%s_%d_init.csv', prm.name, it);
di = readmatrix_ (fullfile (data_dir, 'doe_init', file_grid));
dn = stk_dataframe (di);
zn = stk_dataframe (f (di));
nb_init = size (dn, 1);

Model = stk_model ();
for m = 1:prm.M
    Model(m) = estim_matern (dn, zn(:,m), prm.list_cov, config.lognugget);
end

% candidate grid, common to all maps
xg = stk_sampling_regulargrid (grid_pts, prm.dim_x, prm.BOXx);
sg = stk_sampling_regulargrid (grid_pts, prm.dim_s, prm.BOXs);
cand_set = double (adapt_set (xg, sg));

% fine grids for the overlays: regular in s for display, transformed in
% s (equal-weight quantization of P_S) for the probability computations
xf = stk_sampling_regulargrid (PTS_REF, prm.dim_x, prm.BOXx);
sf = stk_sampling_regulargrid (PTS_REF, prm.dim_s, prm.BOXs);
df_disp = double (adapt_set (xf, sf));
bnd_true = reshape (f (df_disp) <= prm.const(2,1), PTS_REF, PTS_REF)';

df = double (adapt_set (xf, s_trnsf (sf)));
zf = f (df);
true_set = get_true_quantile_set (zf, PTS_REF, PTS_REF, prm.alpha, prm.const);

quantOpt.nbLevels = config.nVar;
quantOpt.useGaussHermite = 1;

stats = [];

for t = 0:(max (compare_steps))

    [xt, st, IS] = integration_grid (Model, dn, zn, prm, config, s_trnsf);
    dt = adapt_set (xt, st);
    mu_n = predict_mean_var (Model, dn, zn, dt, prm);
    estim_current = get_estim (mu_n, prm.const, prm, config);

    if ismember (t, compare_steps)

        % n includes the initial design, as in the paper.
        n = size (dn, 1);
        nb_added = n - nb_init;
        fprintf ('Criterion maps at n = %d (%d added evaluations, %d candidates)...\n', ...
            n, nb_added, size (cand_set, 1));

        tic;
        J = sur_criterion_map (cand_set, Model, dn, zn, dt, IS, ...
            quantOpt, prm, config);
        t_sur = toc;

        tic;
        Psi = eem_criterion_map (cand_set, Model, dn, zn, dt, ...
            estim_current, mu_n, IS, prm, config);
        t_eem = toc;

        fprintf ('  QSI-SUR map: %.1f s, QSI-EEM map: %.1f s\n', t_sur, t_eem);

        rho = spearman_corr (-J, Psi);  % J is minimized, hence the sign
        [~, i_sur] = min (J);
        [~, i_eem] = max (Psi);

        % away from the boundary both criteria are nearly constant and
        % ranks are noise: also correlate where the SUR gain is largest
        act = (J <= pctile (J, 0.1));
        rho_act = spearman_corr (-J(act), Psi(act));

        frac_better_sur = mean (J < J(i_eem));
        frac_better_eem = mean (Psi > Psi(i_sur));
        regret_frac = (J(i_eem) - J(i_sur)) / (max (J) - J(i_sur));

        fprintf (['  Spearman = %.3f (%.3f on the %.0f%% active region)\n' ...
            '  SUR opt (%.2f, %.2f), EEM opt (%.2f, %.2f), ' ...
            'better-than-EEM-opt share = %.3f,\n  SUR regret share = ' ...
            '%.3f\n'], rho, rho_act, 100 * mean (act), ...
            cand_set(i_sur,:), cand_set(i_eem,:), frac_better_sur, ...
            regret_frac);

        stats = [stats; nb_added, n, rho, rho_act, mean(act), ...
            cand_set(i_sur,:), cand_set(i_eem,:), ...
            frac_better_sur, frac_better_eem, regret_frac]; %#ok<AGROW>

        % posterior-mean boundary and plug-in set on the reference grids
        mu_disp = predict_mean_var (Model, dn, zn, df_disp, prm);
        bnd_plug = reshape (mu_disp(:,1) <= prm.const(2,1), PTS_REF, PTS_REF)';
        mu_ref = predict_mean_var (Model, dn, zn, df, prm);
        in_C = (mu_ref(:,1) >= prm.const(1,1)) & (mu_ref(:,1) <= prm.const(2,1));
        plug_set = (mean (reshape (in_C, PTS_REF, PTS_REF), 2) <= prm.alpha)';

        % normalized to [0, 1], with 1 at the best candidate: the SUR
        % criterion is nearly constant, its raw map is unreadable
        J_map = (max (J) - J) / (max (J) - min (J));
        Psi_map = Psi / max (Psi);

        % Keep the archived filenames indexed by added evaluations.
        plot_criterion_map (fig_dir, sprintf ('crit_map_sur_n%03d', nb_added), ...
            J_map, sprintf ('QSI-SUR, n = %d', n), cand_set(i_sur,:), ...
            xg, sg, grid_pts, xf, sf, bnd_true, bnd_plug, true_set, ...
            plug_set, dn, nb_init, prm, opts);

        plot_criterion_map (fig_dir, sprintf ('crit_map_eem_n%03d', nb_added), ...
            Psi_map, sprintf ('QSI-EEM, n = %d', n), cand_set(i_eem,:), ...
            xg, sg, grid_pts, xf, sf, bnd_true, bnd_plug, true_set, ...
            plug_set, dn, nb_init, prm, opts);

    end

    if t == max (compare_steps)
        break
    end

    % one QSI-EEM step
    newpt = get_point_EEM (Model, estim_current, dn, zn, xt, st, ...
        prm.const, IS, prm, config, 0);
    dn = stk_dataframe ([dn; newpt]);
    zn = stk_dataframe ([zn; f(newpt)]);
    for m = 1:prm.M
        Model(m) = estim_matern (dn, zn(:,m), prm.list_cov, config.lognugget);
    end

end

csv_file = fullfile (stats_dir, 'crit_comparison_stats.csv');
writematrix_ (stats, csv_file);
fprintf ('Wrote %s\n', csv_file);
fprintf (['Columns: nb_added, nb_evaluations, spearman, spearman_active, ' ...
    'active_share, x_sur, s_sur, x_eem, s_eem, frac_better_sur, ' ...
    'frac_better_eem, sur_regret_share\n']);

fprintf ('Done.\n');


function [xt, st, IS] = integration_grid (Model, dn, zn, prm, config, s_trnsf)

% importance sampling on X, as in QSI_EEM / QSI_SUR

xt = stk_sampling_randomlhs (config.pts_x, prm.dim_x, prm.BOXx);
st = s_trnsf (stk_sampling_randomlhs (config.pts_s, prm.dim_s, prm.BOXs));

proba_x = proba_tau (Model, dn, zn, xt, st, prm, config);
misclass_x = min (proba_x, 1 - proba_x);

[~, sort_misclass_x] = sort (misclass_x, 'descend');
ranking_x = sort_misclass_x(1);

if sum (misclass_x > 0) <= config.keep_x
    ranking_x = sort_misclass_x(1:config.keep_x);
else
    candidate_set = setdiff (find (misclass_x > 0), ranking_x);
    while size (ranking_x, 1) < config.keep_x
        candidate = randsample (candidate_set, 1, true, ...
            misclass_x(candidate_set));
        candidate_set = setdiff (candidate_set, candidate);
        ranking_x = [ranking_x; candidate]; %#ok<AGROW>
    end
end

xt = xt(ranking_x, :);

IS = 1 ./ misclass_x(ranking_x);
IS(~ isfinite (IS)) = 0;
if sum (IS) > 0
    IS = IS / sum (IS);
else
    IS = ones (size (IS)) / numel (IS);
end

end % function


function J = sur_criterion_map (cand_set, Model, dn, zn, dt, IS, ...
    quantOpt, prm, config)

% same computation as in QSI_SUR.m, with sample paths generated once
% on [dn; cand_set; dt] and rows [dn; pt; dt] extracted per candidate

nb_cand = size (cand_set, 1);
nb_dn = size (dn, 1);
nb_dt = size (dt, 1);

traj_base = zeros (nb_dn + nb_cand + nb_dt, config.nTraj, prm.M);
for m = 1:prm.M
    traj_base(:,:,m) = stk_generate_samplepaths ...
        (Model(m), [dn; cand_set; dt], config.nTraj);
end

J = nan (nb_cand, 1);
traj = zeros (nb_dn + 1 + nb_dt, config.nTraj, prm.M);

for i = 1:nb_cand

    pt = cand_set(i, :);

    for m = 1:prm.M
        traj(1:nb_dn, :, m) = traj_base(1:nb_dn, :, m);
        traj(nb_dn + 1, :, m) = traj_base(nb_dn + i, :, m);
        traj(nb_dn + 2:end, :, m) = traj_base(nb_dn + nb_cand + 1:end, :, m);
    end

    p = stk_predict (Model(1), dn, zn(:,1), pt);
    [var_nodes, weight] = quantization (p.mean, sqrt (p.var), quantOpt);

    bool = cell (config.nVar, prm.M);
    for m = 1:prm.M
        [~, lambda] = stk_predict (Model(m), [dn; pt], [], [dn; pt; dt]);
        lambda_dn = lambda(1:nb_dn, :);
        lambda_pt = lambda(nb_dn + 1, nb_dn + 2:end);

        trajCond = stk_conditioning (lambda_dn, zn(:,m), traj(:,:,m), 1:nb_dn);
        trajCond = trajCond(nb_dn + 2:end, :);

        lambda_pt = reshape (lambda_pt, config.keep_x, config.pts_s);
        tensor_dn = reshape (trajCond, config.keep_x, config.pts_s, config.nTraj);
        traj_ind = reshape (traj(nb_dn + 1, :, m), 1, 1, config.nTraj);

        for k = 1:config.nVar
            bool{k, m} = check_constraints_trajs ...
                (tensor_dn, traj_ind, var_nodes(k), lambda_pt, prm.const(:,m));
        end
    end

    proba = proba_tau_trajs (bool, prm.alpha);
    s = min (proba, 1 - proba);
    J(i) = (IS' * s) * weight(:);

end

end % function


function Psi = eem_criterion_map (cand_set, Model, dn, zn, dt, ...
    estim_current, mu_n, IS, prm, config)

nb_cand = size (cand_set, 1);
Psi = nan (nb_cand, 1);

for i = 1:nb_cand
    Psi(i) = modif_mu (cand_set(i,:), [], Model, dn, zn, dt, ...
        estim_current, mu_n, IS, prm.const, prm, config);
end

end % function


function plot_criterion_map (fig_dir, fig_name, crit, ...
    panel_title, opt_pt, xg, sg, grid_pts, xf, sf, bnd_true, bnd_plug, ...
    true_set, plug_set, dn, nb_init, prm, opts)

fig = paper_new_figure (opts);
ax = paper_new_axes (fig, opts);
set (ax, 'XGrid', 'off', 'YGrid', 'off', 'Layer', 'top');

% a raster image, not pcolor: patch seams show through ps2pdf
image (ax, 'XData', double (xg), 'YData', double (sg), ...
    'CData', reshape (crit, grid_pts, grid_pts)', 'CDataMapping', 'scaled');
colormap (ax, parula);
caxis (ax, [0 1]);
colorbar (ax, 'FontSize', opts.font_size);

contour (ax, double (xf), double (sf)', bnd_true, [1 1], ...
    'LineWidth', 1, 'LineColor', 'black');
contour (ax, double (xf), double (sf)', bnd_plug, [1 1], ...
    'LineWidth', 1, 'LineColor', 'red');

s_lo = min (double (sf));
s_hi = max (double (sf));
plot_set_band (ax, xf, true_set, s_lo + 0.02 * (s_hi - s_lo), 'black');
plot_set_band (ax, xf, plug_set, s_lo + 0.055 * (s_hi - s_lo), 'red');

dn = double (dn);
scatter (ax, dn(1:nb_init, 1), dn(1:nb_init, 2), 12, 'black', 'filled');
if size (dn, 1) > nb_init
    scatter (ax, dn(nb_init + 1:end, 1), dn(nb_init + 1:end, 2), 16, ...
        'red', 'filled', 'MarkerEdgeColor', 'black');
end

plot (ax, opt_pt(1), opt_pt(2), 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', ...
    'LineWidth', 1.5);

xlim (ax, double (prm.BOXx));
ylim (ax, double (prm.BOXs));
xlabel (ax, 'x', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 's', 'FontSize', opts.font_size, 'Color', 'black');
title (ax, panel_title, 'FontSize', opts.font_size, 'FontWeight', 'normal');

paper_save_figure (fig, fig_dir, fig_name);
close (fig);

end % function


function plot_set_band (ax, xf, set_ind, y, color)

sel = nan (numel (set_ind), 1);
sel(set_ind == 1) = y;
line (ax, double (xf(:)), sel, 'LineWidth', 3, 'Color', color);

end % function


function rho = spearman_corr (a, b)

ra = midrank (a);
rb = midrank (b);
ra = ra - mean (ra);
rb = rb - mean (rb);
rho = (ra' * rb) / sqrt ((ra' * ra) * (rb' * rb));

end % function


function q = pctile (x, p)

xs = sort (x(:));
q = xs(max (1, ceil (p * numel (xs))));

end % function


function r = midrank (x)

n = numel (x);
[xs, order] = sort (x(:));
rs = zeros (n, 1);
i = 1;
while i <= n
    j = i;
    while (j < n) && (xs(j+1) == xs(i))
        j = j + 1;
    end
    rs(i:j) = (i + j) / 2;
    i = j + 1;
end
r = zeros (n, 1);
r(order) = rs;

end % function

% Parallel-coordinate figure of the ROTOR37 case (three panels): a
% reference sample of the target set colored by the true violation
% probability, then a sample of the estimated set and the last batch of
% each of the 100 runs, colored by the probability estimated under the
% final GP of the run.  The set samples come from splitting passes on
% the simulator and on the final GPs, rebuilt from the outputs of
% run_small_rotor37.

close all;
fprintf ('ROTOR37 parallel-coordinate plots\n');

script_dir = fileparts (mfilename ('fullpath'));
data_dir   = fullfile (script_dir, 'outputs', 'results');
design_dir = fullfile (data_dir, 'design');
param_dir  = fullfile (data_dir, 'param');
fig_dir    = fullfile (data_dir, 'figures');
if ~ isfolder (fig_dir); mkdir (fig_dir); end

case_name  = 'rotor37';
batch_size = 5;
run_ids    = 1:100;

[problem, config] = paper_get_test_case (case_name);
[prm, f, s_trnsf] = problem ();

opts = struct ( ...
    'fig_width_cm',  13, ...
    'fig_height_cm', 5, ...
    'font_size',     8, ...
    'tick_font_size', 10, ...   % x_n labels
    'nb_ref',        2000, ...   % population of the splitting pass on the simulator
    'nb_meta',       200, ...    % population of the splitting pass on each final GP
    'nb_lines',      500, ...    % lines drawn per plot
    'ratio_range',   [0.6 1.4]); % common color scale for p / alpha

st = stk_sampling_sobol (config.pts_s_deter, prm.dim_s, prm.BOXs);
st = s_trnsf (st);
ns = size (st, 1);


%% Reference sample of the target set, true probability

paper_set_seed (1);
sample_true = set_sample (f, st, prm, config, opts.nb_ref);
p_true = get_proba_deter (eval_product (f, sample_true, st), prm.const, prm, config);
in_true = (p_true <= prm.alpha);
fprintf ('reference: %d points, %d in the target set\n', numel (in_true), sum (in_true));

order = 1:prm.dim_x;
fig = draw_pcp (sample_true(in_true, :), p_true(in_true), ns, order, prm, opts);
paper_save_figure (fig, fig_dir, sprintf ('pcp_%s_true', case_name));
close (fig);


%% Last batches and estimated-set samples, probability under the final GP of each run

x_batch = [];  p_batch = [];
x_meta  = [];  p_meta  = [];

for i = 1:numel (run_ids)
    r = run_ids(i);
    dn = readmatrix_ (fullfile (design_dir, sprintf ('doe_QSI_EEM_smc_%s_%d_%d.csv', ...
        case_name, r, batch_size)));
    zn = f (dn);
    Model = final_model (param_dir, case_name, r, batch_size, prm, config);
    f_pred = @(d) posterior_means (Model, dn, zn, d);

    paper_set_seed (r);
    xm = set_sample (f_pred, st, prm, config, opts.nb_meta);
    xb = dn(end-batch_size+1:end, 1:prm.dim_x);
    p = get_proba_deter (eval_product (f_pred, [xb; xm], st), prm.const, prm, config);

    x_batch = [x_batch; xb];                   %#ok<AGROW>
    p_batch = [p_batch; p(1:batch_size)];      %#ok<AGROW>
    x_meta  = [x_meta; xm];                    %#ok<AGROW>
    p_meta  = [p_meta; p(batch_size+1:end)];   %#ok<AGROW>

    fprintf ('  run %d / %d\n', i, numel (run_ids));
end

fprintf ('last batches: %d points, %d with estimated probability <= alpha\n', ...
    numel (p_batch), sum (p_batch <= prm.alpha));
fig = draw_pcp (x_batch, p_batch, ns, order, prm, opts);
paper_save_figure (fig, fig_dir, sprintf ('pcp_%s_batch', case_name));
close (fig);

in_meta = (p_meta <= prm.alpha);
fprintf ('estimated sets: %d points, %d in the set of their run\n', ...
    numel (in_meta), sum (in_meta));
fig = draw_pcp (x_meta(in_meta, :), p_meta(in_meta), ns, order, prm, opts);
paper_save_figure (fig, fig_dir, sprintf ('pcp_%s_meta', case_name));
close (fig);


%% Local functions

function xt = set_sample (g, st, prm, config, n)

% Population of a deterministic splitting pass on g (same steps as
% deterministic_subset), moved once more inside the final set.

config.pts_x_deter = n;
xt = stk_sampling_maximinlhs (n, prm.dim_x, prm.BOXx);
sigma = config.sigma_init;
thresh = NaN;
first = true;
while true
    if first
        zt = eval_product (g, xt, st);
        [const, thresh] = deterministic_set_const_init (zt, prm, config);
        first = false;
    else
        xt = resample_xt_deter (xt, in_level, config);
        [xt, sigma] = MH_walk_deter (g, xt, st, const, sigma, prm, config);
        zt = eval_product (g, xt, st);
        [const, thresh] = deterministic_set_const (zt, thresh, prm, config);
    end
    in_level = get_estim_deter (zt, const, prm, config);
    if thresh == final_threshold (prm)
        break
    end
    if ~ any (in_level)
        error ('splitting stalled before reaching the target set');
    end
end
xt = resample_xt_deter (xt, in_level, config);
xt = double (MH_walk_deter (g, xt, st, const, sigma, prm, config));

end % function


function Model = final_model (param_dir, case_name, r, batch_size, prm, config)

% Final GP of a run rebuilt from the saved parameters (last row).

Model = stk_model ();
for m = 1:prm.M
    P = readmatrix_ (fullfile (param_dir, sprintf ('param_QSI_EEM_smc_%d_%s_%d_%d.csv', ...
        m, case_name, r, batch_size)));
    idx = readmatrix_ (fullfile (param_dir, sprintf ('cov_QSI_EEM_smc_%d_%s_%d_%d.csv', ...
        m, case_name, r, batch_size)));
    Model(m) = stk_model (prm.list_cov{idx(end)}, prm.dim_x + prm.dim_s);
    Model(m).lognoisevariance = config.lognugget;
    Model(m).param = P(end, :)';
end

end % function


function fig = draw_pcp (X, p, ns, order, prm, opts)

% One line per point, colored by p / alpha: blue ramp inside the set,
% red ramp outside.

n = size (X, 1);
keep = unique (round (linspace (1, n, min (n, opts.nb_lines))));
X = X(keep, order);
v = p(keep) / prm.alpha;

vs = sort (v);
fprintf ('    p / alpha: min %.2f, median %.2f, max %.2f\n', ...
    vs(1), vs(ceil (end / 2)), vs(end));

lo = opts.ratio_range(1);
hi = opts.ratio_range(2);
nb = round (60 * (1 - lo) / (hi - lo));   % colors below alpha
nr = 60 - nb;
wb = linspace (0, 1, nb)';
wr = linspace (0, 1, nr)';
cmap = [(1 - wb) * [0.78 0.86 0.96] + wb * [0.02 0.18 0.40]; ...
        (1 - wr) * [0.98 0.70 0.55] + wr * [0.60 0.05 0.05]];
ci = 1 + round (59 * (min (max (v, lo), hi) - lo) / (hi - lo));

fig = paper_new_figure (opts);
ax = paper_new_axes (fig, opts);
set (ax, 'XGrid', 'off', 'YGrid', 'off', 'FontSize', opts.tick_font_size);

% points near the boundary first, the clearly inside or outside ones on top
[~, idx] = sort (abs (v - 1));
pos = 1:size (X, 2);
for i = idx'
    plot (ax, pos, X(i,:), '-', 'Color', cmap(ci(i),:), 'LineWidth', 0.3);
end

set (ax, 'XLim', [1 pos(end)], 'YLim', [0 1], 'XTick', pos, ...
    'XTickLabel', arrayfun (@(d) sprintf ('x_{%d}', d), order, 'UniformOutput', false));

% color scale drawn as patches: rendered identically by every backend
set (ax, 'Position', [0.06 0.14 0.82 0.83]);
cb = axes ('Parent', fig, 'Position', [0.905 0.14 0.02 0.83], 'FontSize', opts.font_size - 1, ...
    'Box', 'on', 'Layer', 'top', 'XColor', 'black', 'YColor', 'black');
hold (cb, 'on');
edges = linspace (lo, hi, 61);
for k = 1:60
    patch (cb, [0 1 1 0], edges([k k k+1 k+1]), cmap(k,:), 'EdgeColor', 'none');
end
set (cb, 'XLim', [0 1], 'YLim', [lo hi], 'XTick', [], 'YAxisLocation', 'right', ...
    'YTick', lo:(hi - lo) / 4:hi);
ylabel (cb, 'p / \alpha', 'FontSize', opts.font_size, 'Color', 'black');

end % function

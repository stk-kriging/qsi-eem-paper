% ROTOR37 figures of the paper and of the supplement from the campaign
% outputs: scatter of the relative error (same style as the synthetic
% cases), budget convergence, SMC level descent, 2D projections of the
% sampled points, violation probabilities of the last batches and
% per-coordinate profiles.  The parallel-coordinate figure has its own
% script, make_pcp_figures_small_rotor37.  Sections whose input files are
% missing are skipped with a warning (the convergence figure and the set
% samples come from extract_convergence_small_rotor37).

% requires the proprietary ROTOR37 code (see test_cases/rotor37/code)
rotor37_check_available ();

close all;

fprintf('Making ROTOR37 figures and table\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir   = fullfile (script_dir, 'outputs', 'results');
input_dir  = fullfile (data_dir, 'dev');
design_dir = fullfile (data_dir, 'design');
logs_dir   = fullfile (data_dir, 'smc_logs');
param_dir  = fullfile (data_dir, 'param');
fig_dir    = fullfile (data_dir, 'figures');
table_dir  = fullfile (data_dir, 'tables');

for d = {fig_dir, table_dir}
    if ~ isfolder (d{1})
        mkdir (d{1});
    end
end

case_name = 'rotor37';
batch_size = 5;

[problem, config] = paper_get_test_case (case_name);
[prm, f, s_trnsf] = problem ();
pts_init = config.pts_init;

opts = struct ( ...
    'fig_width_cm',  7.5, ...    % width used in the paper
    'fig_height_cm', 5.6, ...
    'font_size',     8, ...
    'nb_bins',       15, ...
    'clip_quantile', 0.99, ...   % upper limit of the relative error axis
    'clip_ratio',    1.25, ...   % clip only if max > clip_ratio * quantile
    'marker',        '*', ...
    'marker_size',   4, ...
    'ref_color',     [0.35 0.35 0.35], ...
    'true_color',    [0.25 0.40 0.75]);

qsi_color = paper_algo_colors ('QSI_EEM');

prefix = sprintf ('QSI_EEM_smc_%s_%d', case_name, batch_size);


%% Scatter of the relative error against the number of batches, and table

d = load_run_data (input_dir, prefix, batch_size);

if isempty (d)
    warning ('scatter: missing dev/size_doe files in %s, skipped.', input_dir);
else
    lim = compute_limits (d, opts);
    fig = draw_panel (d, lim, opts, qsi_color);
    paper_save_figure (fig, fig_dir, sprintf ('scatter_%s_%d', case_name, batch_size));
    close (fig);

    write_table (table_dir, case_name, d);
end


%% Convergence of the deviation with the budget

dev_file = fullfile (input_dir, sprintf ('dev_budgets_%s.txt', prefix));
bud_file = fullfile (input_dir, sprintf ('budgets_%s.txt', prefix));

if ~ (isfile (dev_file) && isfile (bud_file))
    warning (['convergence: missing %s, skipped ' ...
        '(run extract_convergence_small_rotor37 first).'], dev_file);
else
    dev_budgets = readmatrix_ (dev_file);
    budgets = readmatrix_ (bud_file);

    fig = draw_convergence (budgets(:)', dev_budgets, opts, qsi_color);
    paper_save_figure (fig, fig_dir, sprintf ('convergence_%s', case_name));
    close (fig);
end


%% Descent of the SMC tolerance levels

nb_file = fullfile (input_dir, sprintf ('size_doe_%s.txt', prefix));
logs = dir (fullfile (logs_dir, sprintf ('log_%s_*_%d.csv', ...
    strrep (prefix, sprintf ('_%d', batch_size), ''), batch_size)));

if isempty (logs) || ~ isfile (nb_file)
    warning ('levels: missing smc_logs or size_doe file, skipped.');
else
    nb = readmatrix_ (nb_file);
    T = max (nb);
    thresh = nan (numel (logs), T);

    for r = 1:numel (logs)
        L = readmatrix_ (fullfile (logs(r).folder, logs(r).name));
        run_id = sscanf (logs(r).name, ['log_QSI_EEM_smc_' case_name '_%d']);
        for t = 1:min (nb(run_id), T)
            rows = L(:,1) <= t;
            if any (rows)
                thresh(r, t) = L(find (rows, 1, 'last'), 4);
            end
        end
        thresh(r, nb(run_id):T) = prm.tol;  % finished: hold the final level
    end

    fig = draw_levels (thresh, prm.tol, opts, qsi_color);
    paper_save_figure (fig, fig_dir, sprintf ('levels_%s', case_name));
    close (fig);
end


%% Sampled points and set samples: shared inputs of the last three parts

true_file   = fullfile (input_dir, sprintf ('sample_true_%s.txt', case_name));
pooled_file = fullfile (input_dir, sprintf ('sample_meta_pooled_%s.txt', case_name));

designs = dir (fullfile (design_dir, sprintf ('doe_%s_*.csv', prefix(1:end-2))));

sample_true = [];
if isfile (true_file)
    sample_true = readmatrix_ (true_file);
    sample_true = sample_true(:, 1:prm.dim_x);
end

sample_pooled = [];
if isfile (pooled_file)
    sample_pooled = readmatrix_ (pooled_file);
    sample_pooled = sample_pooled(:, 1:prm.dim_x);
end

% activity ranking of the x coordinates, from the fitted GP ranges
params = dir (fullfile (param_dir, sprintf ('param_QSI_EEM_smc_*_%s_*_%d.csv', ...
    case_name, batch_size)));
order = 1:prm.dim_x;
if isempty (params)
    warning ('activity ranking: missing param files, keeping the natural order.');
else
    score = zeros (1, prm.dim_x);
    for r = 1:numel (params)
        P = readmatrix_ (fullfile (params(r).folder, params(r).name));
        score = score + P(end, 2:(1 + prm.dim_x));
    end
    [~, order] = sort (score, 'descend');
end
dims = order(1:3);

% pooled post-initial points of every run, with their batch index
pts = [];
bidx = [];
last_batch = [];
for r = 1:numel (designs)
    D = readmatrix_ (fullfile (designs(r).folder, designs(r).name));
    ns = size (D, 1) - pts_init;
    pts = [pts; D(pts_init+1:end, 1:prm.dim_x)];              %#ok<AGROW>
    bidx = [bidx; ceil((1:ns)' / batch_size)];                %#ok<AGROW>
    last_batch = [last_batch; D(end-batch_size+1:end, 1:prm.dim_x)]; %#ok<AGROW>
end


%% 2D projections on the most active coordinate pairs

pairs = [dims(1) dims(2); dims(1) dims(3); dims(2) dims(3)];

if isempty (sample_true)
    warning (['projections: missing %s, skipped ' ...
        '(run extract_convergence_small_rotor37 first).'], true_file);
elseif isempty (pts)
    warning ('projections: missing design files, skipped.');
else
    % sampled points colored by batch index
    for k = 1:size (pairs, 1)
        fig = draw_projection (sample_true, pts, bidx, pairs(k,:), opts);
        paper_save_figure (fig, fig_dir, sprintf ('proj_%s_%d', case_name, k));
        close (fig);
    end
end


%% Violation probability of the last batches

if isempty (last_batch)
    warning ('violation probabilities: missing design files, skipped.');
else
    st = stk_sampling_sobol (config.pts_s_deter, prm.dim_s, prm.BOXs);
    st = s_trnsf (st);
    zt = eval_product (f, last_batch, st);
    viol = reshape (any (abs (zt) > prm.tol, 2), size (last_batch, 1), []);
    phat = max (mean (viol, 2), 1 / (2 * size (viol, 2)));  % floor at half a count

    fig = draw_phat (phat / prm.alpha, opts, qsi_color);
    paper_save_figure (fig, fig_dir, sprintf ('phat_%s', case_name));
    close (fig);
end


%% Per-coordinate profiles: target set against last batches and estimates

if isempty (sample_true) || isempty (last_batch)
    warning ('profile (last batches): missing inputs, skipped.');
else
    fig = draw_profile (sample_true, last_batch, 'last batches', opts);
    paper_save_figure (fig, fig_dir, sprintf ('profile_%s_last', case_name));
    close (fig);
end

if isempty (sample_true) || isempty (sample_pooled)
    warning ('profile (estimated sets): missing inputs, skipped.');
else
    fig = draw_profile (sample_true, sample_pooled, 'estimated sets', opts);
    paper_save_figure (fig, fig_dir, sprintf ('profile_%s_meta', case_name));
    close (fig);
end


fprintf ('\nDone.\n');


function fig = new_figure (w_cm, h_cm)

fig = figure ('Visible', 'off', 'Color', 'white', 'Renderer', 'painters', ...
    'Units', 'centimeters', 'Position', [2 2 w_cm h_cm], ...
    'PaperUnits', 'centimeters', 'PaperSize', [w_cm h_cm], ...
    'PaperPosition', [0 0 w_cm h_cm]);

% since R2025a figures follow the desktop theme, which prints dark axes
if isprop (fig, 'Theme')
    set (fig, 'Theme', 'light');
end

end % function


function d = load_run_data (input_dir, prefix, batch_size)

dev_file = fullfile (input_dir, ['dev_' prefix '.txt']);
nb_file  = fullfile (input_dir, ['size_doe_' prefix '.txt']);

if ~ (isfile (dev_file) && isfile (nb_file))
    d = [];
    return
end

dev = readmatrix_ (dev_file);
nb  = readmatrix_ (nb_file);

dev = dev(:);
nb  = nb(:);

keep = isfinite (dev) & isfinite (nb);
if ~ all (keep)
    warning ('%s: discarding %d run(s) with a non-finite result.', ...
        prefix, sum (~ keep));
end

d = struct ( ...
    'prefix',     prefix, ...
    'batch_size', batch_size, ...
    'dev',        dev(keep), ...
    'nb',         nb(keep));

end % function


function lim = compute_limits (group, opts)

dev = vertcat (group.dev);
nb  = vertcat (group.nb);

% clip the extreme upper tail, which would flatten the bulk

q = quantile (dev, opts.clip_quantile);

if max (dev) > opts.clip_ratio * q
    lim.clip_level = q;
    lim.clipped = true;
else
    lim.clip_level = max (dev);
    lim.clipped = false;
end

if lim.clip_level <= 0
    lim.clip_level = 1;
end

if lim.clipped
    lim.ylim = [0, 1.09 * lim.clip_level];
else
    lim.ylim = [0, 1.05 * lim.clip_level];
end
lim.yedges = linspace (0, lim.clip_level, opts.nb_bins + 1);

nb_min = min (nb);
nb_max = max (nb);

if nb_max == nb_min
    nb_max = nb_min + 1;
end

if all (nb == round (nb))
    step = max (1, round ((nb_max - nb_min + 1) / opts.nb_bins));
    edges = (nb_min - 0.5):step:(nb_max + 0.5);
    if edges(end) < nb_max + 0.5
        edges(end + 1) = edges(end) + step;
    end
else
    edges = linspace (nb_min, nb_max, opts.nb_bins + 1);
end

lim.xedges = edges;

pad = 0.03 * (edges(end) - edges(1));
lim.xlim = [edges(1) - pad, edges(end) + pad];

end % function


function fig = draw_panel (d, lim, opts, color)

% opaque tint instead of FaceAlpha: painters flattens transparency
face_color = color + (1 - color) * 0.35;

nb  = d.nb;
dev = d.dev;

is_clipped = dev > lim.clip_level;
dev_plot = min (dev, lim.clip_level);

mean_dev = mean (dev);
med_dev  = median (dev);
mean_nb  = mean (nb);
med_nb   = median (nb);

fig = new_figure (opts.fig_width_cm, opts.fig_height_cm);

pos_scatter = [0.350 0.360 0.625 0.610];
pos_hist_y  = [0.040 0.360 0.130 0.610];
pos_hist_x  = [0.350 0.040 0.625 0.130];

ax_x = axes ('Parent', fig, 'Position', pos_hist_x);
histogram (ax_x, nb, lim.xedges, 'Normalization', 'pdf', ...
    'FaceColor', face_color, ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.25);
hold (ax_x, 'on');
ylim_x = get (ax_x, 'YLim');
plot (ax_x, [mean_nb mean_nb], ylim_x, '-', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
plot (ax_x, [med_nb med_nb], ylim_x, '--', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
set (ax_x, 'XLim', lim.xlim, 'YLim', ylim_x, 'YDir', 'reverse', ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none', ...
    'XColor', 'none', 'YColor', 'none');

ax_y = axes ('Parent', fig, 'Position', pos_hist_y);
histogram (ax_y, dev_plot, lim.yedges, 'Normalization', 'pdf', ...
    'Orientation', 'horizontal', ...
    'FaceColor', face_color, ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.25);
hold (ax_y, 'on');
xlim_y = get (ax_y, 'XLim');
plot (ax_y, xlim_y, min (mean_dev, lim.clip_level) * [1 1], '-', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
plot (ax_y, xlim_y, min (med_dev, lim.clip_level) * [1 1], '--', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
set (ax_y, 'XLim', xlim_y, 'YLim', lim.ylim, 'XDir', 'reverse', ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none', ...
    'XColor', 'none', 'YColor', 'none');

ax = axes ('Parent', fig, 'Position', pos_scatter, ...
    'FontSize', opts.font_size, 'Box', 'on', 'Layer', 'top', ...
    'Color', 'white', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

plot (ax, lim.xlim, min (mean_dev, lim.clip_level) * [1 1], '-', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
plot (ax, lim.xlim, min (med_dev, lim.clip_level) * [1 1], '--', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
plot (ax, [mean_nb mean_nb], lim.ylim, '-', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);
plot (ax, [med_nb med_nb], lim.ylim, '--', ...
    'Color', opts.ref_color, 'LineWidth', 0.75);

plot (ax, nb(~ is_clipped), dev_plot(~ is_clipped), opts.marker, ...
    'Color', color, 'MarkerSize', opts.marker_size);

if any (is_clipped)
    plot (ax, nb(is_clipped), dev_plot(is_clipped), '^', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'black', ...
        'MarkerSize', opts.marker_size);
end

set (ax, 'XLim', lim.xlim, 'YLim', lim.ylim);
xlabel (ax, 'Number of batches', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'Relative error', 'FontSize', opts.font_size, 'Color', 'black');

txt = { ...
    'mean / med', ...
    sprintf('err  %.3f / %.3f', mean_dev, med_dev), ...
    sprintf('nb   %.1f / %.1f', mean_nb, med_nb)};

if any (is_clipped)
    txt{end + 1} = sprintf ('%d run(s) > %.3f', ...
        sum (is_clipped), lim.clip_level);
end

annotation (fig, 'textbox', [0.030 0.020 0.305 0.320], 'String', txt, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', opts.font_size - 1, 'Interpreter', 'none', 'Color', 'black', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off', 'Margin', 0);

end % function


function fig = draw_convergence (budgets, dev, opts, color)

qs = quantile (dev, [0.10 0.25 0.50 0.75 0.90], 1);

% log scale: the deviation spans decades over the budget range
floor_val = max (min (dev(dev > 0)) / 2, 1e-4);
qs = max (qs, floor_val);

band_outer = color + (1 - color) * 0.80;
band_inner = color + (1 - color) * 0.55;

fig = new_figure (10, 6.5);
ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Layer', 'top', 'YScale', 'log', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

fill (ax, [budgets fliplr(budgets)], [qs(1,:) fliplr(qs(5,:))], band_outer, ...
    'EdgeColor', 'none');
fill (ax, [budgets fliplr(budgets)], [qs(2,:) fliplr(qs(4,:))], band_inner, ...
    'EdgeColor', 'none');
plot (ax, budgets, qs(3,:), '-o', 'Color', color, 'LineWidth', 1.25, ...
    'MarkerSize', 3.5, 'MarkerFaceColor', color);

set (ax, 'XLim', [budgets(1) budgets(end)]);
xlabel (ax, 'Number of batches', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'Relative error', 'FontSize', opts.font_size, 'Color', 'black');

text (ax, budgets(end), qs(3,end), sprintf ('  %.3f', qs(3,end)), ...
    'FontSize', opts.font_size - 1, 'Color', opts.ref_color, ...
    'VerticalAlignment', 'middle');

end % function


function fig = draw_levels (thresh, tol, opts, color)

T = size (thresh, 2);
qs = quantile (thresh, [0.10 0.50 0.90], 1);

band_color = color + (1 - color) * 0.80;
gray = [0.82 0.82 0.82];

fig = new_figure (10, 6.5);
ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Layer', 'top', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

fill (ax, [1:T, T:-1:1], [qs(1,:), fliplr(qs(3,:))], band_color, ...
    'EdgeColor', 'none');
plot (ax, 1:T, thresh', '-', 'Color', gray, 'LineWidth', 0.25);
plot (ax, 1:T, qs(2,:), '-', 'Color', color, 'LineWidth', 1.25);
plot (ax, [1 T], [tol tol], '--', 'Color', opts.ref_color, 'LineWidth', 0.75);
text (ax, T, tol, sprintf (' tol = %.2f', tol), ...
    'FontSize', opts.font_size - 1, 'Color', opts.ref_color, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

set (ax, 'XLim', [1 T], 'YLim', [tol - 0.005, max(thresh(:)) + 0.01]);
xlabel (ax, 'Number of batches', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'Tolerance level', 'FontSize', opts.font_size, 'Color', 'black');

end % function


function fig = draw_projection (sample_true, pts, bidx, pair, opts)

% at most ~8000 sampled points: enough for the figure, light enough for eps
n = size (pts, 1);
keep = unique (round (linspace (1, n, min (n, 8000))));

nc = max (bidx);
cmap = paper_color_ramp (nc);

fig = new_figure (7.5, 6.2);
ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Layer', 'top', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

plot (ax, sample_true(:,pair(1)), sample_true(:,pair(2)), '.', ...
    'Color', opts.true_color, 'MarkerSize', 4);
scatter (ax, pts(keep,pair(1)), pts(keep,pair(2)), 5, bidx(keep), 'filled');

colormap (ax, cmap);
clim (ax, [1 nc]);
cb = colorbar (ax);
cb.Label.String = 'Batch';
cb.Label.FontSize = opts.font_size;
cb.FontSize = opts.font_size - 1;

axis (ax, [0 1 0 1]);
axis (ax, 'square');
xlabel (ax, sprintf ('x_{%d}', pair(1)), 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, sprintf ('x_{%d}', pair(2)), 'FontSize', opts.font_size, 'Color', 'black');

end % function


function fig = draw_phat (ratio, opts, color)

face_color = color + (1 - color) * 0.35;

fig = new_figure (10, 6.5);
ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Layer', 'top', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

lr = log10 (ratio);
edges = linspace (min (lr) - 0.05, max (max (lr), 0.3) + 0.05, 25);
histogram (ax, lr, edges, 'Normalization', 'pdf', ...
    'FaceColor', face_color, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.25);

yl = get (ax, 'YLim');
plot (ax, [0 0], yl, '--', 'Color', opts.ref_color, 'LineWidth', 0.75);
text (ax, 0, yl(2), sprintf (' %.0f%% in \\Gamma', 100 * mean (ratio <= 1)), ...
    'FontSize', opts.font_size - 1, 'Color', opts.ref_color, ...
    'VerticalAlignment', 'top');

set (ax, 'YLim', yl);
xlabel (ax, 'log_{10} of the violation probability over \alpha', ...
    'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'Density', 'FontSize', opts.font_size, 'Color', 'black');

end % function


function fig = draw_profile (sample_true, sample_other, other_label, opts)

red = [0.75 0.25 0.25];
dim = size (sample_true, 2);

fig = new_figure (9, 7);
ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Layer', 'top', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

for d = 1:dim
    y = dim - d + 1;   % coordinate 1 on top
    qt = quantile (sample_true(:,d), [0.10 0.50 0.90]);
    qo = quantile (sample_other(:,d), [0.10 0.50 0.90]);

    plot (ax, qt([1 3]), [y y] + 0.16, '-', 'Color', opts.true_color, ...
        'LineWidth', 2.5);
    plot (ax, qt(2), y + 0.16, 'o', 'Color', opts.true_color, ...
        'MarkerFaceColor', 'white', 'MarkerSize', 3.5);
    plot (ax, qo([1 3]), [y y] - 0.16, '-', 'Color', red, 'LineWidth', 2.5);
    plot (ax, qo(2), y - 0.16, 'o', 'Color', red, ...
        'MarkerFaceColor', 'white', 'MarkerSize', 3.5);
end

set (ax, 'XLim', [0 1], 'YLim', [0.4 dim + 0.6], 'YTick', 1:dim, ...
    'YTickLabel', arrayfun (@(d) sprintf ('x_{%d}', d), dim:-1:1, ...
    'UniformOutput', false));
xlabel (ax, 'Value', 'FontSize', opts.font_size, 'Color', 'black');
title (ax, sprintf ('target set (blue) vs %s (red), [10-90]%% and median', ...
    other_label), 'FontSize', opts.font_size, 'FontWeight', 'normal', ...
    'Color', 'black');

end % function


function write_table (table_dir, case_name, d)

nb  = d.nb;
pts = nb * d.batch_size;
dev = d.dev;

col_name = sprintf ('QSI-EEM ($r=%d$)', d.batch_size);
nb_row  = fmt_cell (nb, 1, 0);
pts_row = fmt_cell (pts, 1, 0);
err_row = fmt_cell (dev, 3, 3);

tex_file = fullfile (table_dir, sprintf ('table_%s.tex', case_name));

fid = fopen (tex_file, 'wt');
if fid < 0
    error ('Cannot open %s for writing.', tex_file);
end

fprintf (fid, '%% Generated by make_figures_small_rotor37.m\n');
fprintf (fid, '\\begin{table}[ht]\n');
fprintf (fid, '  \\centering\n');
fprintf (fid, '  \\footnotesize\n');
fprintf (fid, '  \\begin{tabular}{|l|c|}\n');
fprintf (fid, '    \\hline\n');
fprintf (fid, '    & %s \\\\\n', col_name);
fprintf (fid, '    \\hline\n');
fprintf (fid, '    No. of batches & %s \\\\\n', nb_row);
fprintf (fid, '    \\hline\n');
fprintf (fid, '    No. of points & %s \\\\\n', pts_row);
fprintf (fid, '    \\hline\n');
fprintf (fid, '    Rel. error & %s \\\\\n', err_row);
fprintf (fid, '    \\hline\n');
fprintf (fid, '  \\end{tabular}\n');
fprintf (fid, ['  \\caption{Summary of QSI-EEM metrics on the ROTOR37 ' ...
    'case. Average (min - max).}\n']);
fprintf (fid, '  \\label{eem:table:rotor37}\n');
fprintf (fid, '\\end{table}\n');

fclose (fid);
fprintf ('Wrote %s\n', tex_file);

fprintf ('\n  %-18s%s\n', sprintf ('rotor37 (n = %d)', numel (dev)), col_name);
fprintf ('  %-18s%s\n', 'No. of batches', strrep (nb_row, '\scriptsize', ''));
fprintf ('  %-18s%s\n', 'No. of points', strrep (pts_row, '\scriptsize', ''));
fprintf ('  %-18s%s\n', 'Rel. error', strrep (err_row, '\scriptsize', ''));

end % function


function s = fmt_cell (v, mean_digits, range_digits)

s = sprintf ('%.*f \\scriptsize(%.*f -- %.*f)', ...
    mean_digits, mean (v), range_digits, min (v), range_digits, max (v));

end % function

% Figures 6.2 to 6.5 and the summary table of Section 6.2, from the data
% written by extract_results_small_synthetic.  One scatter plot of the
% relative error against the number of batches per algorithm, with the two
% marginals as histograms.

close all;

fprintf('Making small synthetic figures and tables\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir  = fullfile (script_dir, 'outputs', 'results');
input_dir = fullfile (data_dir, 'dev');
fig_dir   = fullfile (data_dir, 'figures');
table_dir = fullfile (data_dir, 'tables');

for d = {fig_dir, table_dir}
    if ~ isfolder (d{1})
        mkdir (d{1});
    end
end

case_name = 'piston_smc';

% Order in which the test cases are numbered in Section 6.2
case_order = {'piston_smc', 'trid_smc', 'otl_smc'};
case_index = find (strcmp (case_name, case_order), 1);
case_label = strrep (case_name, '_smc', '');

batch_sizes = [1 2 3];

opts = struct ( ...
    'fig_width_cm',  7.5, ...    % width used in the paper
    'fig_height_cm', 5.6, ...
    'font_size',     8, ...
    'nb_bins',       15, ...
    'clip_quantile', 0.99, ...   % upper limit of the relative error axis
    'clip_ratio',    1.25, ...   % clip only if max > clip_ratio * quantile
    'marker',        '*', ...
    'marker_size',   4, ...
    'ref_color',     [0.35 0.35 0.35]);

qsi_color = paper_algo_colors ('QSI_EEM');
bss_color = paper_algo_colors ('BSS');

fprintf ('\n=== %s (case %d)\n', case_name, case_index);

% QSI-EEM-SMC: one data set per batch size, sharing the same axes

qsi = [];
for b = batch_sizes
    prefix = sprintf ('QSI_EEM_smc_%s_%d', case_name, b);
    d = load_run_data (input_dir, prefix, b);
    if ~ isempty (d)
        qsi = [qsi d];  %#ok<AGROW>
    end
end

if isempty (qsi)
    error ('No QSI-EEM-SMC data found in %s for %s.', input_dir, case_name);
end

lim = compute_limits (qsi, opts);

for k = 1:numel (qsi)
    fig = draw_panel (qsi(k), lim, opts, qsi_color);
    paper_save_figure (fig, fig_dir, ...
        sprintf ('scatter_%s_%d', case_name, qsi(k).batch_size));
    close (fig);
end

% BSS: single data set, on its own scale

prefix = sprintf ('BSS_%s_1', case_name);
bss = load_run_data (input_dir, prefix, 1);

if ~ isempty (bss)
    lim_bss = compute_limits (bss, opts);
    fig = draw_panel (bss, lim_bss, opts, bss_color);
    paper_save_figure (fig, fig_dir, sprintf ('scatter_%s_BSS', case_name));
    close (fig);
else
    warning ('No BSS data found for %s.', case_name);
end

write_table (table_dir, case_label, case_index, [qsi bss]);

fprintf ('\nDone.\n');


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

% extract_precision_small returns NaN when the splitting misses the set
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

% headroom so that the markers on the clipping level stay visible
if lim.clipped
    lim.ylim = [0, 1.09 * lim.clip_level];
else
    lim.ylim = [0, 1.05 * lim.clip_level];
end
lim.yedges = linspace (0, lim.clip_level, opts.nb_bins + 1);

% Number of batches: integer-valued, use integer bin edges.

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

fig = figure ('Visible', 'off', 'Color', 'white', 'Renderer', 'painters', ...
    'Units', 'centimeters', ...
    'Position', [2 2 opts.fig_width_cm opts.fig_height_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', [opts.fig_width_cm opts.fig_height_cm], ...
    'PaperPosition', [0 0 opts.fig_width_cm opts.fig_height_cm]);

% since R2025a figures follow the desktop theme, which prints dark axes
if isprop (fig, 'Theme')
    set (fig, 'Theme', 'light');
end

% histograms outside the tick and axis labels of the scatter plot
pos_scatter = [0.350 0.360 0.625 0.610];
pos_hist_y  = [0.040 0.360 0.130 0.610];
pos_hist_x  = [0.350 0.040 0.625 0.130];

% Marginal distribution of the number of batches (bottom, pointing down)

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

% Marginal distribution of the relative error (left, pointing left)

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

% Scatter plot

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

% Summary of both marginals (solid line: mean, dashed line: median)

txt = { ...
    'mean / med', ...
    sprintf('err  %.3f / %.3f', mean_dev, med_dev), ...
    sprintf('nb   %.1f / %.1f', mean_nb, med_nb)};

if any (is_clipped)
    txt{end + 1} = sprintf ('%d run(s) > %.3f', ...
        sum (is_clipped), lim.clip_level);
end

% lower-left corner, left empty above: text extents are estimated here,
% not measured, so an annotation inside the plot risks being clipped
annotation (fig, 'textbox', [0.030 0.020 0.305 0.320], 'String', txt, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', opts.font_size - 1, 'Interpreter', 'none', 'Color', 'black', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off', 'Margin', 0);

end % function


function write_table (table_dir, case_label, case_index, group)

% one column per element of GROUP, as in Tables 6.1--6.3

col_names = cell (1, numel (group));
for k = 1:numel (group)
    if strncmp (group(k).prefix, 'BSS', 3)
        col_names{k} = 'BSS';
    else
        col_names{k} = sprintf ('QSI-EEM ($r=%d$)', group(k).batch_size);
    end
end

nb_row  = cell (1, numel (group));
pts_row = cell (1, numel (group));
err_row = cell (1, numel (group));

for k = 1:numel (group)
    nb = group(k).nb;
    pts = nb * group(k).batch_size;
    dev = group(k).dev;
    nb_row{k}  = fmt_cell (nb, 1, 0);
    pts_row{k} = fmt_cell (pts, 1, 0);
    err_row{k} = fmt_cell (dev, 3, 3);
end

% LaTeX

col_spec = ['|l|' repmat('c|', 1, numel (group))];
tex_file = fullfile (table_dir, sprintf ('table_%s.tex', case_label));

fid = fopen (tex_file, 'wt');
if fid < 0
    error ('Cannot open %s for writing.', tex_file);
end

fprintf (fid, '%% Generated by make_figures_small_synthetic.m\n');
fprintf (fid, '\\begin{table}[ht]\n');
fprintf (fid, '  \\begin{adjustwidth}{-3cm}{-3cm}\n');
fprintf (fid, '  \\centering\n');
fprintf (fid, '  \\footnotesize\n');
fprintf (fid, '  \\begin{tabular}{%s}\n', col_spec);
fprintf (fid, '    \\hline\n');
fprintf (fid, '    %s \\\\\n', join_cells ('', col_names));
fprintf (fid, '    \\hline\n');
fprintf (fid, '    No. of batches %s \\\\\n', join_cells ('', nb_row));
fprintf (fid, '    \\hline\n');
fprintf (fid, '    No. of added points %s \\\\\n', join_cells ('', pts_row));
fprintf (fid, '    \\hline\n');
fprintf (fid, '    Rel. error %s \\\\\n', join_cells ('', err_row));
fprintf (fid, '    \\hline\n');
fprintf (fid, '  \\end{tabular}\n');
fprintf (fid, ['    \\caption{Summary of QSI-EEM and BSS metrics for ' ...
    'case %d. Average (min - max).}\n'], case_index);
fprintf (fid, '  \\label{eem:table:%s}\n', case_label);
fprintf (fid, '  \\end{adjustwidth}\n');
fprintf (fid, '\\end{table}\n');

fclose (fid);
fprintf ('Wrote %s\n', tex_file);

% Console

nb_txt  = clean_cells (nb_row);
pts_txt = clean_cells (pts_row);
err_txt = clean_cells (err_row);

fprintf ('\n  %-18s',sprintf ('%s (n = %d)', case_label, numel (group(1).dev)));
fprintf ('%-26s', col_names{:});
fprintf ('\n  %-18s','No. of batches');
fprintf ('%-26s', nb_txt{:});
fprintf ('\n  %-18s','No. of points');
fprintf ('%-26s', pts_txt{:});
fprintf ('\n  %-18s','Rel. error');
fprintf ('%-26s', err_txt{:});
fprintf ('\n');

end % function


function s = fmt_cell (v, mean_digits, range_digits)

s = sprintf ('%.*f \\scriptsize(%.*f -- %.*f)', ...
    mean_digits, mean (v), range_digits, min (v), range_digits, max (v));

end % function


function s = join_cells (prefix, cells)

s = prefix;
for k = 1:numel (cells)
    s = [s sprintf('& %s ', cells{k})];  %#ok<AGROW>
end

end % function


function cells = clean_cells (cells)

for k = 1:numel (cells)
    cells{k} = strrep (cells{k}, '\scriptsize', '');
end

end % function

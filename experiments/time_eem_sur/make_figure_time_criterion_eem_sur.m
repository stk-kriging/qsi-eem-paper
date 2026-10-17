% Figure and table of the criterion runtime, from the measurements of
% run_time_criterion_eem_sur.  Reports the cost over one design step: the
% preparation, paid once, plus config.keep_xs candidate evaluations.  Both
% axes logarithmic, so a cost in M^p is a line of slope p.

close all;

fprintf ('Making QSI-SUR / QSI-EEM runtime figure and table\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir  = fullfile (script_dir, 'outputs', 'results');
input_file = fullfile (data_dir, 'time', 'criterion_time.csv');
case_name = 'branin_mod';   % same test case as run_time_criterion_eem_sur
fig_dir   = fullfile (data_dir, 'figures');
table_dir = fullfile (data_dir, 'tables');

for d = {fig_dir, table_dir}
    if ~ isfolder (d{1})
        mkdir (d{1});
    end
end

if ~ isfile (input_file)
    error (['No measurement found in %s.\n' ...
        'Run run_time_criterion_eem_sur first.'], input_file);
end

opts = struct ( ...
    'fig_width_cm',  14.0, ...   % included at 14cm in the paper
    'fig_height_cm', 8.0, ...
    'font_size',     9, ...
    'line_width',    1.5, ...
    'ref_slopes',    [2 3], ...  % reference lines of the log-log plot
    'marker_size',   4, ...
    'sur_color',     paper_algo_colors ('QSI_SUR'), ...
    'eem_color',     paper_algo_colors ('QSI_EEM'), ...
    'all_sizes',     false, ...  % table: every measured grid size
    'table_sizes',   [1000 5000 10000 15000 20000 25000], ...
    'nb_digits',     3);         % table: significant digits of the times

% the preparation is paid once for config.keep_xs evaluations, not for one
[~, config] = paper_get_test_case (case_name, false);
nb_candidates = config.keep_xs;

raw = readmatrix_ (input_file);

grid_size = raw(:,1);
sur_prepare = raw(:,3);
sur_total = raw(:,3) + nb_candidates * raw(:,4);
eem_total = raw(:,5) + nb_candidates * raw(:,6);

sizes = unique (grid_size);
sur_med = nan (size (sizes));
eem_med = nan (size (sizes));
prep_share = nan (size (sizes));

for k = 1:numel (sizes)
    sel = (grid_size == sizes(k));
    sur_med(k) = median (sur_total(sel), 'omitnan');
    eem_med(k) = median (eem_total(sel), 'omitnan');
    prep_share(k) = median (sur_prepare(sel), 'omitnan') / sur_med(k);
end

% a size that could not be measured is dropped, not plotted as zero
sur_ok = isfinite (sur_med) & (sur_med > 0);
eem_ok = isfinite (eem_med) & (eem_med > 0);

if ~ all (sur_ok & eem_ok)
    warning ('%d grid size(s) without a complete measurement.', ...
        sum (~ (sur_ok & eem_ok)));
end

fprintf ('%d grid size(s), from %d to %d, %d candidate(s) per step\n', ...
    numel (sizes), min (sizes), max (sizes), nb_candidates);
fprintf ('  QSI-SUR: %.4g s to %.4g s\n', ...
    min (sur_med(sur_ok)), max (sur_med(sur_ok)));
fprintf ('  QSI-EEM: %.4g s to %.4g s\n', ...
    min (eem_med(eem_ok)), max (eem_med(eem_ok)));

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

ax = axes ('Parent', fig, 'FontSize', opts.font_size, 'Box', 'on', ...
    'Color', 'white', 'XColor', 'black', 'YColor', 'black', ...
    'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], ...
    'GridAlpha', 1);
hold (ax, 'on');

% log-log, ticks in seconds: a cost in M^p is a line of slope p
set (ax, 'XScale', 'log', 'YScale', 'log');

h_sur = plot (ax, sizes(sur_ok), sur_med(sur_ok), ...
    '-o', 'Color', opts.sur_color, 'MarkerFaceColor', opts.sur_color, ...
    'LineWidth', opts.line_width, 'MarkerSize', opts.marker_size);

h_eem = plot (ax, sizes(eem_ok), eem_med(eem_ok), ...
    '-o', 'Color', opts.eem_color, 'MarkerFaceColor', opts.eem_color, ...
    'LineWidth', opts.line_width, 'MarkerSize', opts.marker_size);

x_min = min (sizes);
x_max = max (sizes);
set (ax, 'XLim', [0.9 * x_min, 1.1 * x_max]);

% bound on the data: a steep reference line would flatten both curves
y_lo = min ([sur_med(sur_ok); eem_med(eem_ok)]);
y_hi = max ([sur_med(sur_ok); eem_med(eem_ok)]);

% Reference slopes, anchored on the QSI-SUR curve at the largest grid
h_ref = gobjects (1, numel (opts.ref_slopes));
ref_labels = cell (1, numel (opts.ref_slopes));
x_ref = [x_min, x_max];
for k = 1:numel (opts.ref_slopes)
    q = opts.ref_slopes(k);
    y_ref = max (sur_med(sur_ok)) * (x_ref / x_max) .^ q;
    h_ref(k) = plot (ax, x_ref, y_ref, '--', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.75);
    ref_labels{k} = sprintf ('slope %d', q);
end

set (ax, 'YLim', [0.7 * y_lo, 1.4 * y_hi]);

xlabel (ax, 'grid size', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'time (s)', 'FontSize', opts.font_size, 'Color', 'black');

legend (ax, [h_sur h_eem h_ref], [{'QSI-SUR', 'QSI-EEM'}, ref_labels], ...
    'FontSize', opts.font_size, 'Location', 'northwest');

paper_save_figure (fig, fig_dir, 'time_comp_1');
close (fig);

% table, from the same medians as the figure

if opts.all_sizes
    keep = true (size (sizes));
else
    keep = ismember (sizes, opts.table_sizes(:));
    if ~ any (keep)
        warning ('None of the requested grid sizes was measured: using all.');
        keep = true (size (sizes));
    end
end

keep = keep & sur_ok & eem_ok;

table_sizes = sizes(keep);
table_sur = sur_med(keep);
table_eem = eem_med(keep);
table_share = prep_share(keep);
table_ratio = table_sur ./ table_eem;

nb_repetitions = sum (grid_size == sizes(1));

tex_file = fullfile (table_dir, 'table_time_eem_sur.tex');

fid = fopen (tex_file, 'wt');
if fid < 0
    error ('Cannot open %s for writing.', tex_file);
end

fprintf (fid, '%% Generated by make_figure_time_criterion_eem_sur.m\n');
fprintf (fid, '\\begin{table}[ht]\n');
fprintf (fid, '  \\centering\n');
fprintf (fid, '  \\begin{tabular}{|r|c|c|c|}\n');
fprintf (fid, '    \\hline\n');
fprintf (fid, '    $|\\tX\\times\\tS|$ & QSI-SUR & QSI-EEM & Ratio \\\\\n');
fprintf (fid, '    \\hline\n');

for k = 1:numel (table_sizes)
    fprintf (fid, '    $%d$ & $%s$ & $%s$ & $%s$ \\\\\n', ...
        table_sizes(k), ...
        fmt_time (table_sur(k), opts.nb_digits), ...
        fmt_time (table_eem(k), opts.nb_digits), ...
        fmt_ratio (table_ratio(k)));
end

fprintf (fid, '    \\hline\n');
fprintf (fid, '  \\end{tabular}\n');
fprintf (fid, ['  \\caption{Time (in seconds) spent in the criterion over ' ...
    'one sequential\n    design step, given the size of the approximation ' ...
    'grid $|\\tX\\times\\tS|$:\n    the preparation, paid once, plus %d ' ...
    'candidate evaluations.\n    Median over %d repetitions.}\n'], ...
    nb_candidates, nb_repetitions);
fprintf (fid, '  \\label{eem:table:comp_time}\n');
fprintf (fid, '\\end{table}\n');

fclose (fid);
fprintf ('Wrote %s\n', tex_file);

fprintf ('\n  %10s %14s %14s %10s %14s\n', ...
    'grid size', 'QSI-SUR (s)', 'QSI-EEM (s)', 'ratio', 'SUR prepare');

for k = 1:numel (table_sizes)
    fprintf ('  %10d %14s %14s %10s %13.0f%%\n', ...
        table_sizes(k), ...
        fmt_time (table_sur(k), opts.nb_digits), ...
        fmt_time (table_eem(k), opts.nb_digits), ...
        fmt_ratio (table_ratio(k)), ...
        100 * table_share(k));
end

% Empirical growth, over the whole measured range and over its upper half

report_exponent ('whole range', grid_size, sur_total, eem_total, ...
    min (sizes), max (sizes));

half = sizes(max (1, round (numel (sizes) / 2)));
report_exponent ('upper half ', grid_size, sur_total, eem_total, ...
    half, max (sizes));

fprintf ('\nDone.\n');


function s = fmt_time (t, nb_digits)

% times span several orders: fixed significant digits, not decimals

if ~ isfinite (t)
    s = '--';
    return
end

nb_decimals = max (0, nb_digits - 1 - floor (log10 (abs (t))));
s = sprintf ('%.*f', nb_decimals, t);

end % function


function s = fmt_ratio (r)

if ~ isfinite (r)
    s = '--';
elseif r >= 100
    s = sprintf ('%.0f', r);
else
    s = sprintf ('%.1f', r);
end

end % function


function report_exponent (label, grid_size, sur_total, eem_total, lo, hi)

t_sur_lo = median (sur_total(grid_size == lo), 'omitnan');
t_sur_hi = median (sur_total(grid_size == hi), 'omitnan');
t_eem_lo = median (eem_total(grid_size == lo), 'omitnan');
t_eem_hi = median (eem_total(grid_size == hi), 'omitnan');

p_sur = log (t_sur_hi / t_sur_lo) / log (hi / lo);
p_eem = log (t_eem_hi / t_eem_lo) / log (hi / lo);

fprintf ('\n  empirical growth, %s (%d to %d): QSI-SUR M^%.2f, QSI-EEM M^%.2f\n', ...
    label, lo, hi, p_sur, p_eem);

end % function

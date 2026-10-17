% Figures and Table SM1 of the supplement, sensitivity to the SMC
% parameters, from the data of extract_results_eem_smc_parameters_comparison.
% opts.panels selects the overlay of the three values of a parameter, one
% panel per value, or both.

close all;

fprintf ('Making EEM-SMC parameter comparison figures and table\n');

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

case_name = 'trid_smc';
batch_size = 1;

% name of the example in the supplement
case_labels = struct ('piston_smc', 'piston', ...
    'trid_smc', '$7$-Trid', 'otl_smc', 'OTL');
case_label = case_labels.(case_name);

% same variations as in run_eem_smc_parameters_comparison, with the name
% used in the supplement
variations = { ...
    'rho',            [0.15, 0.35, 0.5], 'Proportion $\rho$'                        ; ...
    'size_stop_test', [100, 300, 500],   'Inducing-set size $|\mathds{O}_n|$'       ; ...
    'stop',           [0.05, 0.1, 0.15], 'Stopping tolerance $\tau$'                ; ...
    'kappa',          [0.05, 0.1, 0.2],  'Relaxation factor $\kappa$'               ; ...
    'MH_steps',       [10, 25, 50],      'Metropolis--Hastings moves'               ; ...
    'MH_accept',      [0.2, 0.35, 0.5],  'Target acceptance rate $A_{\mathrm{target}}$' ; ...
    'pts_x',          [125, 250, 500],   'Particle number $|\tX_1|$'                };

opts = struct ( ...
    'panels',        'overlay', ...  % 'overlay', 'single' or 'both'
    'fig_width_cm',  7.5, ...        % width used in the paper
    'fig_height_cm', 5.6, ...
    'font_size',     8, ...
    'nb_bins',       15, ...
    'clip_quantile', 0.99, ...       % upper limit of the relative error axis
    'clip_ratio',    1.25, ...       % clip only if max > clip_ratio * quantile
    'marker_size',   4, ...
    'markers',       {{'o', '*', 's'}}, ...
    'colors',        paper_color_ramp (3));

fprintf ('\n=== %s\n', case_name);

% read from the configuration, so the table follows the code
[~, config] = paper_get_test_case (case_name, false);

stats = cell (size (variations, 1), 1);

for v = 1:size (variations, 1)

    parameter_name = variations{v, 1};
    parameter_values = variations{v, 2};

    group = [];
    for j = 1:numel (parameter_values)
        d = load_run_data (input_dir, case_name, batch_size, ...
            parameter_name, parameter_values(j), j, config);
        if ~ isempty (d)
            group = [group d];  %#ok<AGROW>
        end
    end

    if isempty (group)
        continue
    end

    fprintf ('  %s\n', parameter_name);

    stats{v} = group;

    % Common axis limits and bin edges for the values of this parameter
    lim = compute_limits (group, opts);

    if any (strcmp (opts.panels, {'overlay', 'both'}))
        fig = draw_panel (group, lim, opts, parameter_name);
        paper_save_figure (fig, fig_dir, sprintf ('scatter_%s_%d__%s', ...
            case_name, batch_size, parameter_name));
        close (fig);
    end

    if any (strcmp (opts.panels, {'single', 'both'}))
        for k = 1:numel (group)
            fig = draw_panel (group(k), lim, opts, parameter_name);
            paper_save_figure (fig, fig_dir, sprintf ('scatter_%s_%d__%s_%s', ...
                case_name, batch_size, parameter_name, group(k).value_str));
            close (fig);
        end
    end

end

if all (cellfun (@isempty, stats))
    error ('No data found in %s for %s.', input_dir, case_name);
end

write_table (table_dir, case_name, case_label, variations, stats);

fprintf ('\nDone.\n');


function d = load_run_data (input_dir, case_name, batch_size, ...
    parameter_name, value, rank, config)

% the value is part of the file name, as num2str writes it
value_str = num2str (value);

suffix = sprintf ('QSI_EEM_smc_%s_%d__%s_%s.csv', ...
    case_name, batch_size, parameter_name, value_str);

dev_file = fullfile (input_dir, ['dev_' suffix]);
nb_file  = fullfile (input_dir, ['size_doe_' suffix]);

if ~ (isfile (dev_file) && isfile (nb_file))
    d = [];
    return
end

dev = readmatrix_ (dev_file);
nb  = readmatrix_ (nb_file);

dev = dev(:);
nb  = nb(:);

% NaN when the splitting misses the quantile set: not plottable
keep = isfinite (dev) & isfinite (nb);
if ~ all (keep)
    warning ('%s = %s: discarding %d run(s) with a non-finite result.', ...
        parameter_name, value_str, sum (~ keep));
end

d = struct ( ...
    'parameter',  parameter_name, ...
    'value',      value, ...
    'value_str',  value_str, ...
    'rank',       rank, ...
    'is_default', value == config.(parameter_name), ...
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


function fig = draw_panel (group, lim, opts, parameter_name)

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
pos_scatter = [0.350 0.360 0.625 0.550];
pos_hist_y  = [0.040 0.360 0.130 0.550];
pos_hist_x  = [0.350 0.040 0.625 0.130];

ax_x = axes ('Parent', fig, 'Position', pos_hist_x);
hold (ax_x, 'on');

ax_y = axes ('Parent', fig, 'Position', pos_hist_y);
hold (ax_y, 'on');

ax = axes ('Parent', fig, 'Position', pos_scatter, ...
    'FontSize', opts.font_size, 'Box', 'on', 'Layer', 'top', ...
    'Color', 'white', 'XColor', 'black', 'YColor', 'black');
hold (ax, 'on');

nb_clipped = 0;
h_scatter = gobjects (1, numel (group));
labels = cell (1, numel (group));

for k = 1:numel (group)

    % colour and marker follow the rank, so panels match the overlay
    idx = 1 + mod (group(k).rank - 1, size (opts.colors, 1));
    color = opts.colors(idx, :);
    marker = opts.markers{1 + mod (idx - 1, numel (opts.markers))};

    nb = group(k).nb;
    dev = group(k).dev;

    is_clipped = dev > lim.clip_level;
    dev_plot = min (dev, lim.clip_level);
    nb_clipped = nb_clipped + sum (is_clipped);

    % marginals as outlines, so that they can be superimposed
    histogram (ax_x, nb, lim.xedges, 'Normalization', 'pdf', ...
        'DisplayStyle', 'stairs', 'EdgeColor', color, 'LineWidth', 0.75);
    histogram (ax_y, dev_plot, lim.yedges, 'Normalization', 'pdf', ...
        'Orientation', 'horizontal', ...
        'DisplayStyle', 'stairs', 'EdgeColor', color, 'LineWidth', 0.75);

    h_scatter(k) = plot (ax, nb(~ is_clipped), dev_plot(~ is_clipped), ...
        marker, 'Color', color, 'MarkerSize', opts.marker_size, ...
        'LineWidth', 0.75);

    if any (is_clipped)
        plot (ax, nb(is_clipped), dev_plot(is_clipped), '^', ...
            'MarkerFaceColor', color, 'MarkerEdgeColor', 'black', ...
            'MarkerSize', opts.marker_size);
    end

    % mean of the two marginals
    plot (ax, mean (nb), min (mean (dev), lim.clip_level), 'd', ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'black', ...
        'MarkerSize', opts.marker_size + 2, 'LineWidth', 0.5);

    % value only: the parameter name overflowed the box, which then vanished
    labels{k} = group(k).value_str;
    if group(k).is_default
        labels{k} = [labels{k} ' (default)'];
    end

end

set (ax, 'XLim', lim.xlim, 'YLim', lim.ylim);
xlabel (ax, 'Number of batches', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'Relative error', 'FontSize', opts.font_size, 'Color', 'black');
title (ax, strrep (parameter_name, '_', '\_'), 'FontSize', opts.font_size, ...
    'Color', 'black', 'Interpreter', 'tex');

xlim_y = get (ax_y, 'XLim');
set (ax_y, 'XLim', xlim_y, 'YLim', lim.ylim, 'XDir', 'reverse', ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none', ...
    'XColor', 'none', 'YColor', 'none');

ylim_x = get (ax_x, 'YLim');
set (ax_x, 'XLim', lim.xlim, 'YLim', ylim_x, 'YDir', 'reverse', ...
    'XTick', [], 'YTick', [], 'Box', 'off', 'Color', 'none', ...
    'XColor', 'none', 'YColor', 'none');

% lower-left corner, left empty above: text extents are estimated here,
% not measured, so a legend inside the plot risks being clipped
lgd = legend (ax, h_scatter, labels, 'FontSize', opts.font_size - 2, ...
    'Interpreter', 'tex', 'Box', 'off', 'Color', 'none');
set (lgd, 'Units', 'normalized');
set (lgd, 'Position', [0.030 0.090 0.290 0.230]);

if nb_clipped > 0
    annotation (fig, 'textbox', [0.030 0.010 0.290 0.070], ...
        'String', sprintf ('%d run(s) > %.3f', nb_clipped, lim.clip_level), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', opts.font_size - 2, 'Interpreter', 'none', ...
        'Color', 'black', 'EdgeColor', 'none', ...
        'FitBoxToText', 'off', 'Margin', 0);
end

end % function


function write_table (table_dir, case_name, case_label, variations, stats)

short_name = strrep (case_name, '_smc', '');

tex_file = fullfile (table_dir, ...
    sprintf ('table_smc_sensitivity_%s.tex', short_name));

fid = fopen (tex_file, 'wt');
if fid < 0
    error ('Cannot open %s for writing.', tex_file);
end

fprintf (fid, '%% Generated by make_figures_eem_smc_parameters_comparison.m\n');
fprintf (fid, '\\begin{table}[h]\n');
fprintf (fid, '  \\centering\n');
fprintf (fid, ['  \\caption{Relative error and number of additional ' ...
    'evaluations after\n    the initial design in the SMC parameter ' ...
    'sensitivity analysis, on the\n    %s example.  Values are reported ' ...
    'as mean (standard deviation).\n    Default parameter values are ' ...
    'shown in bold.}\n'], case_label);
fprintf (fid, '  \\label{eem:tab:smc_sensitivity_%s}\n', short_name);
fprintf (fid, '  \\begin{tabular}{lccc}\n');
fprintf (fid, '    \\hline\n');
fprintf (fid, ['    Parameter & Value & Relative error & ' ...
    'Additional evaluations \\\\\n']);
fprintf (fid, '    \\hline\n');

fprintf ('\n  %-30s %8s %18s %16s\n', ...
    case_name, 'Value', 'Relative error', 'Add. evals');

for v = 1:size (variations, 1)

    group = stats{v};
    if isempty (group)
        continue
    end

    for k = 1:numel (group)

        dev = group(k).dev;
        nb = group(k).nb;

        value_tex = fmt_value (group(k));
        if group(k).is_default
            value_tex = sprintf ('\\mathbf{%s}', value_tex);
        end

        dev_tex = sprintf ('%.3f\\ (%.3f)', mean (dev), std (dev));
        nb_tex  = sprintf ('%.1f\\ (%.1f)', mean (nb), std (nb));

        if k == 1
            name_tex = variations{v, 3};
        else
            name_tex = '  ';
        end

        fprintf (fid, '    %s & $%s$ & $%s$ & $%s$ \\\\\n', ...
            name_tex, value_tex, dev_tex, nb_tex);

        if k == 1
            name_txt = variations{v, 1};
        else
            name_txt = '';
        end

        fprintf ('  %-30s %8s %8.3f (%.3f) %7.1f (%.1f)\n', ...
            name_txt, fmt_value (group(k)), ...
            mean (dev), std (dev), mean (nb), std (nb));

    end

    fprintf (fid, '    \\hline\n');

end

fprintf (fid, '  \\end{tabular}\n');
fprintf (fid, '\\end{table}\n');

fclose (fid);
fprintf ('\nWrote %s\n', tex_file);

end % function


function s = fmt_value (d)

% the supplement reports the relaxation factor as 1 + kappa
value = d.value;
if strcmp (d.parameter, 'kappa')
    value = 1 + value;
end

if value == round (value)
    s = sprintf ('%d', value);
else
    s = sprintf ('%.2f', value);
end

end % function

% Quadrature size figures of the supplement, from the deviations written
% by extract_results_quadrature_comparison.  The curves lie closer than
% the width of the strokes that draw them, so a companion _rel figure
% holds the same curves relative to the finest quadrature.

close all;

fprintf ('Making quadrature size comparison figures\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir  = fullfile (script_dir, 'outputs', 'results');
input_dir = fullfile (data_dir, 'deviations');
fig_dir   = fullfile (data_dir, 'figures');

if ~ isfolder (fig_dir)
    mkdir (fig_dir);
end

case_name = 'double_camel';
algo_name = 'QSI_EEM';
parameter = 'nVar';
parameter_values = [3, 7, 11, 15, 20];
nb_runs = 100;

opts = struct ( ...
    'fig_width_cm',  7.0, ...   % \toto = 70mm in the paper
    'fig_height_cm', 5.25, ...
    'font_size',     8, ...
    'line_width',    1.0, ...
    'rel_floor',     0.05);   % reference below this fraction of its own
                              % maximum: relative difference not reported

% quantile order, prefix of the figure name
summaries = { ...
    0.50, 'metric' ; ...
    0.75, 'metric_75' ; ...
    0.95, 'metric_95' };

[problem, config] = paper_get_test_case (case_name);
prm = problem ();

steps = 0:config.axT:config.T;

fprintf ('\n=== %s: %d step(s), %s in %s\n', ...
    prm.name, numel (steps), parameter, mat2str (parameter_values));

nb_values = numel (parameter_values);
dev = cell (nb_values, 1);

for k = 1:nb_values
    dev{k} = load_deviations (input_dir, algo_name, prm.name, ...
        parameter, parameter_values(k), nb_runs);
    if isempty (dev{k})
        fprintf ('  %s = %-3d no data, skipped\n', parameter, parameter_values(k));
    else
        fprintf ('  %s = %-3d %d run(s)\n', parameter, parameter_values(k), ...
            size (dev{k}, 1));
    end
end

if all (cellfun (@isempty, dev))
    error ('No data found in %s for %s.', input_dir, prm.name);
end

colors = paper_color_ramp (nb_values);

% Reference of the relative figures: the finest quadrature that has data
ref = find (~ cellfun (@isempty, dev), 1, 'last');

for s = 1:size (summaries, 1)

    level = summaries{s,1};

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);

    h = gobjects (1, nb_values);
    labels = {};
    kept = [];

    for k = 1:nb_values
        if isempty (dev{k})
            continue
        end
        y = quantile (dev{k}, level, 1);
        h(k) = plot (ax, steps, y, '-', 'Color', colors(k,:), ...
            'LineWidth', opts.line_width);
        labels{end+1} = sprintf ('%s = %d', parameter, parameter_values(k));  %#ok<SAGROW>
        kept(end+1) = k;                                                      %#ok<SAGROW>
    end

    set (ax, 'XLim', [steps(1), steps(end)]);
    xlabel (ax, 'steps', 'FontSize', opts.font_size, 'Color', 'black');
    ylabel (ax, 'prop. misclass', 'FontSize', opts.font_size, 'Color', 'black');

    legend (ax, h(kept), labels, 'FontSize', opts.font_size - 1, ...
        'Location', 'northeast', 'Interpreter', 'none');

    paper_save_figure (fig, fig_dir, ...
        sprintf ('%s_%s_%s', summaries{s,2}, prm.name, parameter));
    close (fig);

    % same curves, relative to the finest quadrature

    y_ref = quantile (dev{ref}, level, 1);

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);

    plot (ax, [steps(1), steps(end)], [0 0], '-', ...
        'Color', [0.35 0.35 0.35], 'LineWidth', 0.75);

    h = gobjects (1, nb_values);
    labels = {};
    kept = [];

    for k = 1:nb_values
        if isempty (dev{k}) || (k == ref)
            continue
        end
        y = quantile (dev{k}, level, 1);
        % the ratio explodes where the reference approaches zero
        rel = 100 * (y - y_ref) ./ y_ref;
        rel(y_ref <= opts.rel_floor * max (y_ref)) = NaN;
        h(k) = plot (ax, steps, rel, '-', 'Color', colors(k,:), ...
            'LineWidth', opts.line_width);
        labels{end+1} = sprintf ('%s = %d', parameter, parameter_values(k));  %#ok<SAGROW>
        kept(end+1) = k;                                                      %#ok<SAGROW>
    end

    set (ax, 'XLim', [steps(1), steps(end)]);
    xlabel (ax, 'steps', 'FontSize', opts.font_size, 'Color', 'black');
    ylabel (ax, sprintf ('rel. diff. to %s = %d (%%)', ...
        parameter, parameter_values(ref)), ...
        'FontSize', opts.font_size, 'Color', 'black');

    % 'best': the curves change corner on each test case
    legend (ax, h(kept), labels, 'FontSize', opts.font_size - 1, ...
        'Location', 'best', 'Interpreter', 'none');

    paper_save_figure (fig, fig_dir, ...
        sprintf ('%s_%s_%s_rel', summaries{s,2}, prm.name, parameter));
    close (fig);

end

fprintf ('\nDone.\n');


function dev = load_deviations (input_dir, algo_name, case_name, ...
    parameter, value, nb_runs)

dev = [];

for it = 1:nb_runs
    f = fullfile (input_dir, sprintf ('dev_%s_%s_%d__%s_%d.csv', ...
        algo_name, case_name, it, parameter, value));
    if ~ isfile (f)
        continue
    end
    row = readmatrix_ (f);
    dev = [dev; row(:)'];  %#ok<AGROW>
end

end % function

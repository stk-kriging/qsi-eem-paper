% Figures 4.2 and 4.3: one curve per algorithm, from the deviations
% written by extract_results_qsi_comparison.

close all;

fprintf ('Making QSI comparison figures\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir  = fullfile (script_dir, 'outputs', 'results');
input_dir = fullfile (data_dir, 'deviations');
fig_dir   = fullfile (data_dir, 'figures');

if ~ isfolder (fig_dir)
    mkdir (fig_dir);
end

case_name = 'branin_mod';
nb_runs = 100;

% name in the data files, legend, colour (shared by all the figures)
algos = { ...
    'random',            'random',    paper_algo_colors('random')            ; ...
    'Ranjan',            'Ranjan',    paper_algo_colors('Ranjan')            ; ...
    'misclassification', 'misclass.', paper_algo_colors('misclassification') ; ...
    'joint_m',           'joint-SUR', paper_algo_colors('joint_m')           ; ...
    'QSI_m',             'QSI-SUR',   paper_algo_colors('QSI_m')             ; ...
    'QSI_EEM',           'QSI-EEM',   paper_algo_colors('QSI_EEM')           };

opts = struct ( ...
    'fig_width_cm',  7.0, ...   % \toto = 70mm in the paper
    'fig_height_cm', 5.25, ...
    'font_size',     8, ...
    'line_width',    1.5);

% quantile order, prefix of the figure name
summaries = { ...
    0.50, 'metric' ; ...
    0.75, 'metric_75' ; ...
    0.95, 'metric_95' };

[problem, config] = paper_get_test_case (case_name);
prm = problem ();

steps = 0:config.axT:config.T;

fprintf ('\n=== %s: %d step(s), up to %d run(s) per algorithm\n', ...
    prm.name, numel (steps), nb_runs);

dev = cell (size (algos, 1), 1);
for a = 1:size (algos, 1)
    dev{a} = load_deviations (input_dir, algos{a,1}, prm.name, nb_runs);
    if isempty (dev{a})
        fprintf ('  %-18s no data, skipped\n', algos{a,1});
    else
        fprintf ('  %-18s %d run(s)\n', algos{a,1}, size (dev{a}, 1));
    end
end

if all (cellfun (@isempty, dev))
    error ('No data found in %s for %s.', input_dir, prm.name);
end

for s = 1:size (summaries, 1)

    level = summaries{s,1};

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);

    h = gobjects (1, size (algos, 1));
    labels = {};
    kept = [];

    for a = 1:size (algos, 1)
        if isempty (dev{a})
            continue
        end
        y = quantile (dev{a}, level, 1);
        h(a) = plot (ax, steps, y, '-', 'Color', algos{a,3}, ...
            'LineWidth', opts.line_width);
        labels{end+1} = algos{a,2};  %#ok<SAGROW>
        kept(end+1) = a;             %#ok<SAGROW>
    end

    set (ax, 'XLim', [steps(1), steps(end)]);
    xlabel (ax, 'steps', 'FontSize', opts.font_size, 'Color', 'black');
    ylabel (ax, 'prop. misclass', 'FontSize', opts.font_size, 'Color', 'black');

    % 'best': with six algorithms a fixed corner covers the curves
    legend (ax, h(kept), labels, 'FontSize', opts.font_size - 1, ...
        'Location', 'best', 'Interpreter', 'none');

    paper_save_figure (fig, fig_dir, sprintf ('%s_%s', summaries{s,2}, prm.name));
    close (fig);

end

fprintf ('\nDone.\n');


function dev = load_deviations (input_dir, algo_name, case_name, nb_runs)

dev = [];

for it = 1:nb_runs
    f = fullfile (input_dir, ...
        sprintf ('dev_%s_%s_%d.csv', algo_name, case_name, it));
    if ~ isfile (f)
        continue
    end
    row = readmatrix_ (f);
    dev = [dev; row(:)'];  %#ok<AGROW>
end

end % function

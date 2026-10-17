% Per-algorithm figures of the supplement, median with quantile bands and
% individual sample paths, from the deviations written by
% extract_results_qsi_comparison.  Grey everywhere and a shared vertical
% axis: the six panels are read side by side.

close all;

fprintf ('Making QSI comparison supplementary figures\n');

script_dir = fileparts (mfilename ('fullpath'));

data_dir  = fullfile (script_dir, 'outputs', 'results');
input_dir = fullfile (data_dir, 'deviations');
fig_dir   = fullfile (data_dir, 'figures');

if ~ isfolder (fig_dir)
    mkdir (fig_dir);
end

case_name = 'branin_mod';
nb_runs = 100;

% name in the data files, title, colour (shared by all the figures)
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
    'line_width',    1.5, ...
    'traj_width',    0.25, ...
    'outer_band',    [0.05 0.95], ...
    'inner_band',    [0.25 0.75], ...
    'grey_base',     [0.15 0.15 0.15], ...  % bands and trajectories
    'outer_tint',    0.82, ...   % 1 is white, 0 is grey_base
    'inner_tint',    0.60, ...
    'traj_tints',    [0.90 0.20], ...   % pale to saturated, never darker
    'traj_clip',     0.995);       % upper limit of the trajectory panels

[problem, config] = paper_get_test_case (case_name);
prm = problem ();

steps = 0:config.axT:config.T;

fprintf ('\n=== %s: %d step(s), up to %d run(s) per algorithm\n', ...
    prm.name, numel (steps), nb_runs);

% Load every algorithm first: the vertical axis of each family is shared
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

band_top = max (cellfun (@(d) top_quantile (d, opts.outer_band(2)), dev));
traj_top = max (cellfun (@(d) top_quantile_all (d, opts.traj_clip), dev));
metric_ylim = [0, 1.05 * band_top];
trajs_ylim  = [0, 1.02 * traj_top];

for a = 1:size (algos, 1)

    if isempty (dev{a})
        continue
    end

    color = algos{a,3};
    fig_suffix = sprintf ('%s_%s', prm.name, algos{a,1});

    % median and quantiles

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);

    draw_band (ax, steps, dev{a}, opts.outer_band, ...
        tint (opts.grey_base, opts.outer_tint));
    draw_band (ax, steps, dev{a}, opts.inner_band, ...
        tint (opts.grey_base, opts.inner_tint));

    plot (ax, steps, quantile (dev{a}, 0.5, 1), '-', ...
        'Color', color, 'LineWidth', opts.line_width);

    finish_axes (ax, steps, opts, algos{a,2}, metric_ylim);
    paper_save_figure (fig, fig_dir, ['metric_' fig_suffix]);
    close (fig);

    % individual sample paths

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);

    % sorted by final value: pale where the run converged, dark where not
    [~, order] = sort (dev{a}(:,end), 'ascend');
    traj = dev{a}(order, :);

    set (ax, 'ColorOrder', tint_ramp (opts.grey_base, opts.traj_tints, size (traj, 1)));
    plot (ax, steps, traj', '-', 'LineWidth', opts.traj_width);

    % white halo, so that the median reads through the bundle
    med = quantile (dev{a}, 0.5, 1);
    plot (ax, steps, med, '-', 'Color', 'white', 'LineWidth', 3);
    plot (ax, steps, med, '-', 'Color', color, 'LineWidth', 1.75);

    finish_axes (ax, steps, opts, algos{a,2}, trajs_ylim);
    set (ax, 'Layer', 'bottom');   % grid behind the pale trajectories
    paper_save_figure (fig, fig_dir, ['trajs_' fig_suffix]);
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


function c = tint (color, amount)

c = color + (1 - color) * amount;

end % function


function c = tint_ramp (color, range, n)

% always toward white: the ramp starts from a neutral grey
a = linspace (range(1), range(2), n)';
c = color + a * (1 - color);

end % function


function m = top_quantile (dev, p)

if isempty (dev)
    m = -inf;
else
    m = max (quantile (dev, p, 1));
end

end % function


function m = top_quantile_all (dev, p)

% a few early excursions would otherwise squeeze the bundle downward

if isempty (dev)
    m = -inf;
else
    m = quantile (dev(:), p);
end

end % function


function draw_band (ax, steps, dev, levels, face_color)

lo = quantile (dev, levels(1), 1);
hi = quantile (dev, levels(2), 1);

fill (ax, [steps, fliplr(steps)], [lo, fliplr(hi)], face_color, ...
    'EdgeColor', 'none');

end % function


function finish_axes (ax, steps, opts, title_str, ylim_)

set (ax, 'XLim', [steps(1), steps(end)], 'YLim', ylim_, 'Layer', 'top');
xlabel (ax, 'steps', 'FontSize', opts.font_size, 'Color', 'black');
ylabel (ax, 'prop. misclass', 'FontSize', opts.font_size, 'Color', 'black');
title (ax, title_str, 'FontSize', opts.font_size, 'Color', 'black', ...
    'Interpreter', 'none');

end % function

% Demo of QSI-EEM & related algorithms
%
% Runs QSI-EEM, and optionally the comparison algorithms of contrib-qsi, on
% the modified Branin example from Figure 1. Plots the points selected by
% QSI-EEM and the resulting estimate of the quantile set, and, when several
% algorithms are run, compares the designs they produce.
%
% The test case must satisfy dim_x + dim_s = 2 so that the designs can be
% plotted. The algorithms and the number of evaluations are set by the
% variables algorithms and nb_evaluations below.
%
% Outputs are written to `outputs/`.

close all;

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);
root_dir = fileparts(paper_code_dir);

data_dir = fullfile(script_dir,'outputs');
paper_prepare_data_dir(data_dir);

% debug_mode = true: shortened configuration (see paper_debug_config)
debug_mode = false;

% seeds the RNG and names the initial design file
run_id = 1;

% any two-dimensional case known to paper_get_test_case: 'branin_mod'
case_name = 'branin_mod';

% Algorithms to run: QSI-EEM alone by default. Add any of 'QSI-SUR',
% 'Joint-SUR', 'Ranjan', 'Misclassification', 'Random' (contrib-qsi) to
% compare the designs they produce.
algorithms = {'QSI-EEM'};

% Number of sequential evaluations performed by each algorithm after the
% initial design (one point per step; batches of size r > 1 are only
% available in QSI_EEM_SMC).
nb_evaluations = 10;

[problem, config] = paper_get_test_case (case_name, debug_mode);
config.T = nb_evaluations;
config.axT = 1;
[prm, ~, ~] = problem();
if (prm.dim_x + prm.dim_s) ~= 2
  error(['qsieempaper_demo: the test case must satisfy dim_x + dim_s = 2 ' ...
    'for plotting (got dim_x = %d, dim_s = %d)'], prm.dim_x, prm.dim_s);
end
paper_generate_initial_designs (problem, config, run_id, data_dir);
n_init = config.pts_init;

% known algorithms: label, prefix of the output files, runner
all_labels = {'QSI-EEM','QSI-SUR','Joint-SUR','Ranjan', ...
  'Misclassification','Random'};
all_prefixes = {'QSI_EEM', sprintf('QSI_%s',config.critName), ...
  sprintf('joint_%s',config.critName), 'Ranjan', 'misclassification', 'random'};
all_runners = { ...
  @() QSI_EEM (problem, config, run_id, data_dir, 0), ...
  @() QSI_SUR (problem, config, run_id, data_dir, 0), ...
  @() joint_SUR (problem, config, run_id, data_dir), ...
  @() Ranjan (problem, config, run_id, data_dir), ...
  @() misclassification (problem, config, run_id, data_dir), ...
  @() random (problem, config, run_id, data_dir)};

fprintf('Running Figure 1 algorithm demo on %s\n',case_name);
fprintf('  Algorithms: %s\n',strjoin(algorithms,', '));
fprintf('  Evaluations per algorithm: %d\n',nb_evaluations);

% check the list of algorithms before running anything
indices = zeros(size(algorithms));
for k = 1:numel(algorithms)
  idx = find(strcmp(all_labels,algorithms{k}));
  if isempty(idx)
    error('qsieempaper_demo: unknown algorithm ''%s''',algorithms{k});
  end
  indices(k) = idx;
end
if ~any(indices == 1)
  error('qsieempaper_demo: the list of algorithms must include ''QSI-EEM''');
end
prefixes = all_prefixes(indices);

for k = 1:numel(algorithms)
  time_file = fullfile(data_dir,'results','time', ...
    sprintf('time_%s_%s_%d.csv',prefixes{k},prm.name,run_id));
  paper_set_seed(run_id);  % same seed for every algorithm, as in the experiments
  paper_run_timed_algorithm(algorithms{k},all_runners{indices(k)},time_file);
end

fig_timer = tic();
fig_dir = fullfile(data_dir,'results','figures');
png_file = fullfile(fig_dir,'qsi_eem_uncertainty_demo.png');
plot_succeeded = true;
try
    qsieempaper_demo_plot
    fprintf ('  %-20s done in %.2f s\n','Figure generation', toc (fig_timer));
catch plot_error
    plot_succeeded = false;
    png_file = '[PLOT FAILED]';
    fprintf ('  Figure generation    failed in %.2f s\n', toc (fig_timer));
    fprintf ('  Completed algorithm outputs are still available in %s\n', data_dir);
    fprintf ('  %s\n', plot_error.message);
end

manifest = cell (7, 1);
manifest{1} = sprintf ('debug_mode=%d', debug_mode);
manifest{2} = sprintf ('problem=%s', case_name);
manifest{3} = sprintf ('run_id=%d', run_id);
manifest{4} = sprintf ('evaluations=%d', nb_evaluations);
manifest{5} = sprintf ('initial_design_size=%d', n_init);
manifest{6} = sprintf ('algorithms=%s', strjoin (algorithms, ','));
manifest{7} = sprintf ('figure=%s', png_file);
paper_write_manifest (data_dir, manifest);

fprintf('\nCompleted Figure 1 algorithm demo.\n');
fprintf('  Initial design points: %d\n',n_init);
fprintf('  Evaluations per algorithm: %d\n',nb_evaluations);
if plot_succeeded
  fprintf('  Figure: %s\n',png_file);
  if numel(algorithms) > 1
    fprintf('  Design comparison: %s\n', ...
      fullfile(fig_dir,'qsi_eem_design_comparison.png'));
  end
else
  fprintf('  Figure: not generated\n');
end

cd(root_dir);

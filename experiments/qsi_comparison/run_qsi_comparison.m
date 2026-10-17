% QSI criteria comparison: QSI-EEM against QSI-SUR, joint-SUR, Ranjan,
% misclassification and random sampling, for the "Results on moderately
% difficult problems" section of the paper. The default is a short
% Branin modified debug run.
%
% Outputs are written to `outputs/`.

close all;

fprintf('Launching QSI criteria comparison\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
case_names = {'branin_mod', 'double_camel', 'hart4'};
algorithm_names = {'QSI_EEM', 'QSI_SUR', 'joint_SUR', 'Ranjan', ...
  'misclassification', 'random'};
POOL = 50;

if debug_mode
  run_ids = 1;
  case_names = {'branin_mod'};
  POOL = 2;
end

delete (gcp ('nocreate'));
[clust, clustClean] = get_cluster (); %create temporary files to avoid bug when using several pools.
parpool (clust, POOL);
pctRunOnAll run startup.m
pctRunOnAll warning ('off','all');

fprintf('%d case(s), %d run(s), %d algorithm(s) {%s}, %d workers\n', ...
  numel(case_names), numel(run_ids), numel(algorithm_names), ...
  strjoin(algorithm_names, ', '), POOL);

for p = 1:numel(case_names)
  [problem, config] = paper_get_test_case (case_names{p}, debug_mode);

  fprintf('\n=== %s (%d/%d): %d runs x %d algorithms\n', ...
    case_names{p}, p, numel(case_names), numel(run_ids), numel(algorithm_names));

  paper_generate_initial_designs (problem, config, run_ids, data_dir);

  t_case = tic;

  for a = 1:numel(algorithm_names)
    algo = algorithm_names{a};
    fprintf('  %s | %d runs\n', algo, numel(run_ids));
    t_blk = tic;

    parfor i = 1:numel(run_ids)
      run_id = run_ids(i);
      paper_set_seed(run_id)

      switch algo
        case 'QSI_EEM',           QSI_EEM (problem, config, run_id, data_dir, 0);
        case 'QSI_SUR',           QSI_SUR (problem, config, run_id, data_dir, 0);
        case 'joint_SUR',         joint_SUR (problem, config, run_id, data_dir);
        case 'Ranjan',            Ranjan (problem, config, run_id, data_dir);
        case 'misclassification', misclassification (problem, config, run_id, data_dir);
        case 'random',            random (problem, config, run_id, data_dir);
      end
    end

    fprintf('  %s | done in %.1f min\n', algo, toc(t_blk)/60);
  end

  fprintf('=== %s done in %.1f min\n', case_names{p}, toc(t_case)/60);
end

manifest = cell (4, 1);
manifest{1} = sprintf ('debug_mode=%d', debug_mode);
manifest{2} = sprintf ('run_ids=%s', mat2str (run_ids));
manifest{3} = sprintf ('problems=%s', strjoin (case_names, ','));
manifest{4} = sprintf ('algorithms=%s', strjoin (algorithm_names, ','));
paper_write_manifest (data_dir, manifest);

cd(paper_code_dir);


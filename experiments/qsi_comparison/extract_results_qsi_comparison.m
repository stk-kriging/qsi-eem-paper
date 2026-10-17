% Deviations of the QSI-criteria comparison, through EXTRACT_DEVIATION
% on a fixed product grid, for the figures of the "Results on moderately
% difficult problems" section of the paper.

fprintf('Extracting QSI criteria comparison numerical results\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
results_dir = fullfile(data_dir, 'results');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
case_names = {'branin_mod', 'hart4', 'double_camel'};
algos_names = {"QSI_EEM", "QSI_m", "joint_m","Ranjan", "misclassification", "random"};
POOL = 50; %number of cores to use

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

fprintf('%d case(s), %d run(s), %d algorithm(s), %d workers\n', ...
    numel(case_names), numel(run_ids), numel(algos_names), POOL);

for p = 1:numel(case_names)
    [problem, config] = paper_get_test_case (case_names{p}, debug_mode);

    fprintf('\n=== %s (%d/%d): %d runs x %d algorithms\n', ...
        case_names{p}, p, numel(case_names), numel(run_ids), numel(algos_names));

    for j = 1:numel(algos_names)
        fprintf('  %s | %d runs\n', algos_names{j}, numel(run_ids));
        mean_plugin = strcmp (algos_names{j}, 'QSI_EEM');
        t_blk = tic;

        parfor i = 1:numel(run_ids)
            run_id = run_ids(i);
            paper_set_seed(run_id)
          
            extract_deviation (problem, config, ...
                algos_names{j}, run_id, data_dir, mean_plugin);
            
        end

        fprintf('  %s | done in %.1f min\n', algos_names{j}, toc(t_blk)/60);
    end
end

cd(paper_code_dir);


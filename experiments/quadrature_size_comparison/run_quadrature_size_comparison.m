% Effect of the quadrature size (nVar) on QSI-EEM: campaigns for the
% quadrature figures of the supplement.

fprintf('Launching quadrature size comparison\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
results_dir = fullfile(data_dir, 'results');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
parameter_name = 'nVar';
parameter_values = [3, 7, 11, 15, 20]; % values for config.nVar
case_names = {'branin_mod', 'hart4', 'double_camel'};
POOL = 50; %number of cores to use

if debug_mode
    run_ids = 1;
    parameter_values = [3, 7];
    case_names = {'branin_mod'};
    POOL = 2;
end

delete (gcp ('nocreate'));
[clust, clustClean] = get_cluster (); %create temporary files to avoid bug when using several pools.
parpool (clust, POOL);
pctRunOnAll run startup.m
pctRunOnAll warning ('off','all');

fprintf('%d case(s), %d run(s), %s in %s, %d workers\n', ...
    numel(case_names), numel(run_ids), parameter_name, ...
    mat2str(parameter_values), POOL);

for p = 1:numel(case_names)
    [problem, config] = paper_get_test_case (case_names{p}, debug_mode);

    fprintf('\n=== %s (%d/%d): %d runs x %d %s values\n', ...
        case_names{p}, p, numel(case_names), numel(run_ids), ...
        numel(parameter_values), parameter_name);

    paper_generate_initial_designs (problem, config, run_ids, data_dir);

    for j = 1:numel(parameter_values)
        fprintf('  %s = %g | %d runs\n', ...
            parameter_name, parameter_values(j), numel(run_ids));
        t_blk = tic;

        config.(parameter_name) = parameter_values(j);

        parfor i = 1:numel(run_ids)
            run_id = run_ids(i);
            paper_set_seed(run_id);
            QSI_EEM (problem, config, run_id, data_dir, 0);
        end

        paper_rename_files_with_param (results_dir, '__', ...
            char (parameter_name), num2str (parameter_values(j)));

        fprintf('  %s = %g | done in %.1f min\n', ...
            parameter_name, parameter_values(j), toc(t_blk)/60);
    end
end

cd (paper_code_dir);


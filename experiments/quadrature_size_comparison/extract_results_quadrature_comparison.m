% Deviations of the quadrature-size campaigns (supplement figures),
% through EXTRACT_DEVIATION on a fixed product grid.

fprintf('Extracting quadrature size comparison numerical results\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
results_dir = fullfile(data_dir, 'results');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
parameter = "nVar";
parameter_values = [3, 7, 11, 15, 20];
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
    numel(case_names), numel(run_ids), parameter, ...
    mat2str(parameter_values), POOL);

for p = 1:numel(case_names)
    [problem, config] = paper_get_test_case (case_names{p}, debug_mode);

    fprintf('\n=== %s (%d/%d): %d runs x %d %s values\n', ...
        case_names{p}, p, numel(case_names), numel(run_ids), ...
        numel(parameter_values), parameter);

    for j = 1:numel(parameter_values)
        fprintf('  %s = %g | %d runs\n', ...
            parameter, parameter_values(j), numel(run_ids));
        t_blk = tic;

        parfor i = 1:numel(run_ids)
            run_id = run_ids(i);
            paper_set_seed(run_id);
            
            paper_extract_deviation_variable_params (                      ...
                problem, config, "QSI_EEM", parameter,               ...
                num2str (parameter_values(j)), run_id, data_dir, true);            
        end

        fprintf('  %s = %g | done in %.1f min\n', ...
            parameter, parameter_values(j), toc(t_blk)/60);
    end
end

cd(paper_code_dir);

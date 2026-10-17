% Sensitivity of QSI-EEM-SMC to its parameters: one campaign per varied
% parameter value, for the supplement figures and Table SM1.

close all;

fprintf('Launching QSI-EEM-SMC parameter comparison\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
results_dir = fullfile(data_dir, 'results');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:50;
case_names = {'piston_smc', 'otl_smc', 'trid_smc'};
batch_size = 1;
POOL = feature('numcores') - 1; %number of cores to use

variations = { ...
    'rho',            [0.15, 0.35, 0.5] ; ...  % width of the SMC levels
    'size_stop_test', [100, 300, 500]   ; ...  % test set of the stopping criterion
    'stop',           [0.05, 0.1, 0.15] ; ...  % tolerance (stop and stop_final)
    'kappa',          [0.05, 0.1, 0.2]  ; ...  % relaxation of the safe estimator
    'MH_steps',       [10, 25, 50]      ; ...  % quality of the particle moves
    'MH_accept',      [0.2, 0.35, 0.5]  ; ...  % target acceptance rate of the walk
    'pts_x',          [125, 250, 500]   };     % number of particles

if debug_mode
    run_ids = 1;
    case_names = {'trid_smc', 'piston_smc'}; % two cases: exercises the reset below
    POOL = 1; %number of cores to use

    % rho and kappa are not overridden by paper_debug_config, so both
    % blocks keep a baseline
    variations = { ...
        'rho',            [0.15, 0.35, 0.5] ; ...
        'kappa',          [0.05, 0.1, 0.2]  };
end

delete (gcp ('nocreate'));
[clust, clustClean] = get_cluster (); %create temporary files to avoid bug when using several pools.
pool = parpool (clust, POOL);
par_opts = parforOptions (pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
pctRunOnAll run startup.m
pctRunOnAll warning ('off','all');

nb_blocks = sum(cellfun(@numel, variations(:,2)));

fprintf('%d case(s), %d run(s), batch size %d, %d parameter(s), %d blocks per case, %d workers\n', ...
    numel(case_names), numel(run_ids), batch_size, size(variations, 1), ...
    nb_blocks, POOL);

t_all = tic;

for p = 1:numel(case_names)

    [problem, config_base] = paper_get_test_case (case_names{p}, debug_mode);
    [prm, ~, ~] = problem();

    nb_baseline = sum(cellfun(@(n, vals) any(vals == config_base.(n)), ...
        variations(:,1), variations(:,2)));
    nb_computed = nb_blocks - max(nb_baseline - 1, 0);

    fprintf('\n=== %s (%d/%d): dim_x %d, dim_s %d, const %s, alpha %g\n', ...
        case_names{p}, p, numel(case_names), prm.dim_x, prm.dim_s, ...
        mat2str(double(prm.const), 5), prm.alpha);
    fprintf('    %d blocks (%d computed), %d runs each\n', ...
        nb_blocks, nb_computed, numel(run_ids));

    paper_generate_initial_designs (problem, config_base, run_ids, data_dir);

    t_case = tic;

    baseline_tag = {}; % per case: kept, it would copy the previous case

    for v = 1:size(variations, 1)
        parameter_name = variations{v, 1};
        parameter_values = variations{v, 2};

        fprintf('\n  === %s (%d/%d) in %s\n', parameter_name, v, ...
            size(variations, 1), mat2str(parameter_values));

        for j = 1:numel(parameter_values)
            value = parameter_values(j);
            is_baseline = (value == config_base.(parameter_name));

            if is_baseline && ~isempty(baseline_tag)
                fprintf('    %s = %g | baseline, copied from %s = %s\n', ...
                    parameter_name, value, baseline_tag{1}, baseline_tag{2});

                paper_copy_files_with_param (results_dir, ...
                    baseline_tag{1}, baseline_tag{2}, ...
                    parameter_name, num2str (value));
                continue
            end

            config = config_base;
            config.(parameter_name) = value;
            if strcmp (parameter_name, 'stop')
                config.stop_final = value;
            end

            fprintf('    %s = %g | %d runs\n', parameter_name, value, numel(run_ids));
            t_blk = tic;

            parfor (i = 1:numel(run_ids), par_opts)
                run_id = run_ids(i);
                paper_set_seed(run_id);

                QSI_EEM_SMC (problem, config, run_id, batch_size, data_dir, 1);
            end

            paper_rename_files_with_param (results_dir, '__', ...
                parameter_name, num2str (value));

            if is_baseline
                baseline_tag = {parameter_name, num2str(value)};
            end

            fprintf('    %s = %g | done in %.1f min\n', ...
                parameter_name, value, toc(t_blk)/60);
        end
    end

    fprintf('\n=== %s done in %.1f min\n', case_names{p}, toc(t_case)/60);

end

fprintf('\n=== %d case(s) done in %.1f min\n', numel(case_names), toc(t_all)/60);

cd (paper_code_dir);

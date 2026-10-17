% Deviations of the SMC-parameters campaigns (supplement figures and
% Table SM1), through PAPER_EXTRACT_PRECISION_SMALL_VARIABLE_PARAMS.

close all;

fprintf('Extracting QSI-EEM-SMC parameter comparison numerical results\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:50;
case_names = {'piston_smc', 'otl_smc', 'trid_smc'};
batch_size = 1;
POOL = feature('numcores') - 1; %number of cores to use

variations = { ...
    'rho',            [0.15, 0.35, 0.5] ; ...
    'size_stop_test', [100, 300, 500]   ; ...
    'stop',           [0.05, 0.1, 0.15] ; ...
    'kappa',          [0.05, 0.1, 0.2]  ; ...
    'MH_steps',       [10, 25, 50]      ; ...
    'MH_accept',      [0.2, 0.35, 0.5]  ; ...
    'pts_x',          [125, 250, 500]   };

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

    [problem, config] = paper_get_test_case (case_names{p}, debug_mode);
    [prm, ~, ~] = problem();

    % extract_precision_small itself reads rho, MH_steps and MH_accept: every block is
    % measured with the default config, not the varied one.

    nb_baseline = sum(cellfun(@(n, vals) any(vals == config.(n)), ...
        variations(:,1), variations(:,2)));
    nb_computed = nb_blocks - max(nb_baseline - 1, 0);

    fprintf('\n=== %s (%d/%d): const %s, pts_x_deter %d\n', ...
        case_names{p}, p, numel(case_names), ...
        mat2str(double(prm.const), 5), config.pts_x_deter);
    fprintf('    %d blocks (%d computed), %d runs each\n', ...
        nb_blocks, nb_computed, numel(run_ids));

    t_case = tic;

    baseline_res = {}; % per case: kept, it would report the previous case
    baseline_tag = {};

    for v = 1:size(variations, 1)
        parameter_name = variations{v, 1};
        parameter_values = variations{v, 2};

        fprintf('\n  === %s (%d/%d) in %s\n', parameter_name, v, ...
            size(variations, 1), mat2str(parameter_values));

        for j = 1:numel(parameter_values)
            value = parameter_values(j);
            value_str = num2str (value);
            is_baseline = (value == config.(parameter_name));

            % same DoE, config and seeds: the deviation is already known
            if is_baseline && ~isempty(baseline_res)
                fprintf('    %s = %g: baseline, reused from %s = %s\n', ...
                    parameter_name, value, baseline_tag{1}, baseline_tag{2});

                deviation = baseline_res{1};
                nb_batches = baseline_res{2};
            else
                fprintf('    %s = %g: %d runs\n', parameter_name, value, numel(run_ids));
                t_blk = tic;

                deviation = zeros(numel(run_ids), 1);
                nb_batches = zeros(numel(run_ids), 1);

                parfor (i = 1:numel(run_ids), par_opts)
                    run_id = run_ids(i);
                    paper_set_seed(run_id)

                    [deviation(i), nb_batches(i)] = paper_extract_precision_small_variable_params ( ...
                        problem, config, "QSI_EEM_smc", parameter_name, value_str, ...
                        run_id, batch_size, data_dir);
                end

                fprintf('    %s = %g: done in %.1f min\n', ...
                    parameter_name, value, toc(t_blk)/60);

                if is_baseline
                    baseline_res = {deviation, nb_batches};
                    baseline_tag = {parameter_name, value_str};
                end
            end

            file_deviation = sprintf("dev_QSI_EEM_smc_%s_%d__%s_%s.csv", ...
                prm.name, batch_size, parameter_name, value_str);
            writematrix_(deviation, fullfile(data_dir, "results/dev", file_deviation));
            file_nb_batches = sprintf("size_doe_QSI_EEM_smc_%s_%d__%s_%s.csv", ...
                prm.name, batch_size, parameter_name, value_str);
            writematrix_(nb_batches, fullfile(data_dir, "results/dev", file_nb_batches));

            fprintf('    %s = %g: deviation %.3g +/- %.2g, %.1f evaluations\n', ...
                parameter_name, value, mean(deviation, 'omitnan'), ...
                std(deviation, 'omitnan'), mean(nb_batches, 'omitnan'));
        end
    end

    fprintf('\n=== %s done in %.1f min\n', case_names{p}, toc(t_case)/60);

end

fprintf('\n=== %d case(s) done in %.1f min\n', numel(case_names), toc(t_all)/60);

cd (paper_code_dir);

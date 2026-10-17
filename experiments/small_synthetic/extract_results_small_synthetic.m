% Deviations of the synthetic benchmark, by subset simulation
% (EXTRACT_PRECISION_SMALL), for the scatter figures of the "Synthetic
% functions" section of the paper.

close all;

fprintf('Extracting small synthetic numerical results\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
case_names = {'piston_smc', 'otl_smc', 'trid_smc'};
batch_sizes = [1 2 3];
POOL = feature('numcores') - 1; %number of cores to use

if debug_mode
    run_ids = 1;
    case_names = {'trid_smc'};
    batch_sizes = 1;
    POOL = 1; %number of cores to use
end

delete (gcp ('nocreate'));
[clust, clustClean] = get_cluster (); %create temporary files to avoid bug when using several pools.
pool = parpool (clust, POOL);
par_opts = parforOptions (pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
pctRunOnAll run startup.m

fprintf('%d case(s), %d run(s), batch sizes %s, %d workers\n', ...
    numel(case_names), numel(run_ids), mat2str(batch_sizes), POOL);

for p = 1:numel(case_names)
    [problem, config] = paper_get_test_case (case_names{p}, debug_mode);
    [prm, ~, ~] = problem();

    fprintf('\n=== %s (%d/%d): const %s, pts_x_deter %d\n', ...
        case_names{p}, p, numel(case_names), ...
        mat2str(double(prm.const), 5), config.pts_x_deter);

    % QSI-EEM-SMC
    for b = batch_sizes
        fprintf('    QSI-EEM-SMC batch %d: %d runs\n', b, numel(run_ids));

        deviation = zeros(numel(run_ids), 1);
        nb_batches = zeros(numel(run_ids), 1);

        parfor (i = 1:numel(run_ids), par_opts)
            run_id = run_ids(i);
            paper_set_seed(run_id)

            [deviation(i), nb_batches(i)] = ...
                extract_precision_small(problem, config, "QSI_EEM_smc", run_id, b, data_dir);
        end

        file_deviation = sprintf("dev_QSI_EEM_smc_%s_%d", prm.name, b);
        writematrix_(deviation, fullfile(data_dir, "results/dev", file_deviation));
        file_nb_batches = sprintf("size_doe_QSI_EEM_smc_%s_%d", prm.name, b);
        writematrix_(nb_batches, fullfile(data_dir, "results/dev", file_nb_batches));
    end

end

cd (paper_code_dir);


% Benchmark of the synthetic cases (piston, OTL, Trid): 100 runs of
% QSI-EEM-SMC (batch sizes 1, 2, 3), feeding the scatter figures of the
% "Synthetic functions" section of the paper.

close all;

fprintf('Launching small synthetic benchmark (QSI-EEM-SMC)\n');

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
    [prm, ~, ~] = problem ();

    fprintf('\n=== %s (%d/%d): dim_x %d, dim_s %d, const %s, alpha %g\n', ...
        case_names{p}, p, numel(case_names), prm.dim_x, prm.dim_s, ...
        mat2str(double(prm.const), 5), prm.alpha);
    fprintf('    %d runs, batch sizes %s, stop %g / %g, size_stop_test %d\n', ...
        numel(run_ids), mat2str(batch_sizes), config.stop, config.stop_final, ...
        config.size_stop_test);

    paper_generate_initial_designs (problem, config, run_ids, data_dir);

    t_case = tic;

    for b = 1:numel(batch_sizes)
        bs = batch_sizes(b);
        fprintf('  QSI-EEM-SMC | batch %d | %d runs\n', bs, numel(run_ids));
        t_blk = tic;

        parfor (i = 1:numel(run_ids), par_opts)
            run_id = run_ids(i);
            paper_set_seed (run_id);
            QSI_EEM_SMC (problem, config, run_id, bs, data_dir,1);
        end

        fprintf('  QSI-EEM-SMC | batch %d | done in %.1f min\n', bs, toc(t_blk)/60);
    end

    fprintf('=== %s done in %.1f min\n', case_names{p}, toc(t_case)/60);

end

cd (paper_code_dir);


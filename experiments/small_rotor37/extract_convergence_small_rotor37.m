% Deviations of the ROTOR37 campaign at intermediate budgets
% (EXTRACT_PRECISION_BUDGETS), for the convergence figure, plus reference
% samples of the target and estimated sets for the projection figures.

% requires the proprietary ROTOR37 code (see test_cases/rotor37/code)
rotor37_check_available ();

close all;

fprintf('Extracting ROTOR37 convergence results\n');

script_dir = fileparts(mfilename('fullpath'));
paper_code_dir = fileparts(script_dir);

data_dir = fullfile(script_dir, 'outputs');
paper_prepare_data_dir(data_dir);

debug_mode = false;

run_ids = 1:100;
batch_size = 5;
budgets = [10 15 20 25 30 40 60 90];
POOL = feature('numcores') - 1; %number of cores to use

if debug_mode
    run_ids = 1;
    budgets = [2 4];
    POOL = 1; %number of cores to use
end

delete (gcp ('nocreate'));
[clust, clustClean] = get_cluster (); %create temporary files to avoid bug when using several pools.
pool = parpool (clust, POOL);
par_opts = parforOptions (pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
pctRunOnAll run startup.m

fprintf('%d run(s), batch size %d, budgets %s, %d workers\n', ...
    numel(run_ids), batch_size, mat2str(budgets), POOL);

[problem, config] = paper_get_test_case ('rotor37', debug_mode);
[prm, ~, ~] = problem();

fprintf('\n=== %s: tol %g, alpha %g, pts_x_deter %d\n', ...
    prm.name, prm.tol, prm.alpha, config.pts_x_deter);

nb_runs = numel(run_ids);
deviations = zeros(nb_runs, numel(budgets));
xs_true = cell(nb_runs, 1);
xs_pred = cell(nb_runs, 1);

parfor (i = 1:nb_runs, par_opts)
    run_id = run_ids(i);
    paper_set_seed(run_id)

    [deviations(i, :), ~, xs_true{i}, xs_pred{i}] = extract_precision_budgets ( ...
        problem, config, "QSI_EEM_smc", run_id, batch_size, budgets, data_dir);
end

writematrix_(deviations, fullfile(data_dir, "results/dev", ...
    sprintf("dev_budgets_QSI_EEM_smc_%s_%d", prm.name, batch_size)));
writematrix_(budgets(:), fullfile(data_dir, "results/dev", ...
    sprintf("budgets_QSI_EEM_smc_%s_%d", prm.name, batch_size)));

% reference samples of the target and estimated sets, from the first run
writematrix_(double(xs_true{1}), ...
    fullfile(data_dir, "results/dev", sprintf("sample_true_%s", prm.name)));
writematrix_(double(xs_pred{1}), ...
    fullfile(data_dir, "results/dev", sprintf("sample_meta_%s", prm.name)));

% pooled sample of the estimated sets, ~100 points from every run
pooled = [];
for i = 1:nb_runs
    xp = double(xs_pred{i});
    if ~isempty(xp)
        idx = unique(round(linspace(1, size(xp, 1), min(size(xp, 1), 100))));
        pooled = [pooled; xp(idx, :)]; %#ok<AGROW>
    end
end
writematrix_(pooled, ...
    fullfile(data_dir, "results/dev", sprintf("sample_meta_pooled_%s", prm.name)));

for k = 1:numel(budgets)
    fprintf('    budget %3d batches: deviation %.3g +/- %.2g\n', budgets(k), ...
        mean(deviations(:,k), 'omitnan'), std(deviations(:,k), 'omitnan'));
end

cd (paper_code_dir);

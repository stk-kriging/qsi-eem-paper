% Called from qsieempaper_demo, whose workspace provides data_dir,
% case_name, run_id and debug_mode, and optionally the cell arrays
% algorithms (labels) and prefixes (output file prefixes) of the algorithms
% that were run: when several are present, a second figure compares their
% designs.

[problem, config] = paper_get_test_case (case_name, debug_mode);
[prm, f, s_trnsf] = problem();

init_file = fullfile(data_dir,'doe_init', ...
  sprintf('doe_init_%s_%d_init.csv',prm.name,run_id));
initial_design = readmatrix_(init_file);
n_init = size(initial_design,1);

design_file = fullfile(data_dir,'results','design', ...
  sprintf('doe_QSI_EEM_%s_%d.csv',prm.name,run_id));
design = readmatrix_(design_file);
max_step = size(design,1) - n_init;
if max_step < 0
  error('QSI-EEM design has fewer rows than the initial design.');
end

param_tables = cell(prm.M,1);
cov_tables = cell(prm.M,1);
max_saved_step = max_step;
for m = 1:prm.M
  param_file = fullfile(data_dir,'results','param', ...
    sprintf('param_QSI_EEM_%d_%s_%d.csv',m,prm.name,run_id));
  cov_file = fullfile(data_dir,'results','param', ...
    sprintf('cov_QSI_EEM_%d_%s_%d.csv',m,prm.name,run_id));
  param_tables{m} = readmatrix_(param_file);
  cov_tables{m} = readmatrix_(cov_file);
  max_saved_step = min(max_saved_step, size(param_tables{m},1) - 1);
  max_saved_step = min(max_saved_step, size(cov_tables{m},1) - 1);
end

steps = unique([0 max_saved_step]);
n_steps = numel(steps);

pts_x = 180;
pts_s = 180;
xf = stk_sampling_regulargrid(pts_x,prm.dim_x,prm.BOXx);
x = double(xf);
x = x(:);

sf_regular = stk_sampling_regulargrid(pts_s,prm.dim_s,prm.BOXs);
s_regular = double(sf_regular);
s_regular = s_regular(:);
df_regular = adapt_set(xf,sf_regular);
zf_regular = f(df_regular);
z_grid = reshape(zf_regular,pts_x,pts_s)';
c_grid = double(reshape((zf_regular <= prm.const(2,1)),pts_x,pts_s)');

sf = stk_sampling_sobol(pts_s,prm.dim_s,prm.BOXs);
sf = s_trnsf(sf);
sf = sort(sf);
df = adapt_set(xf,sf);
[true_set,true_proba] = get_true_quantile_set( ...
  f(df),pts_x,pts_s,prm.alpha,prm.const);

uncert_config = config;
uncert_config.pts_x = pts_x;
uncert_config.pts_s = pts_s;
uncert_config.nTraj = max(config.nTraj,50);

fig = figure('visible','off');
clf(fig);
set(fig,'paperpositionmode','auto');
set(fig,'position',[100,100,520*n_steps,700]);

for k = 1:n_steps
  step = steps(k);
  row_index = step + 1;
  dn = design(1:(n_init + step),:);
  zn = f(dn);

  Model = stk_model();
  for m = 1:prm.M
    cov_index = cov_tables{m}(row_index,:);
    cov = prm.list_cov{cov_index};
    Model(m) = stk_model(cov,prm.dim_x + prm.dim_s);
    Model(m).param = param_tables{m}(row_index,:);
  end

  [estimated_set, estimated_proba] = get_expected_quantile_set( ...
    Model,df,pts_x,pts_s,dn,zn,prm.const,prm.alpha);
  set_error = lebesgue_deviation(true_set,estimated_set);

  paper_set_seed(9000 + step);
  class_proba = proba_tau(Model,dn,zn,xf,sf,prm,uncert_config);
  class_uncertainty = min(class_proba,1 - class_proba);

  subplot(3,n_steps,k);
  p = pcolor(x,s_regular,z_grid);
  set(p,'EdgeColor','none');
  hold on;
  contour(x,s_regular,c_grid,[1 1],'black','LineWidth',1.5);
  plot(initial_design(:,1),initial_design(:,2),'.', ...
    'Color',[0.25 0.25 0.25],'MarkerSize',12);
  if step > 0
    added = design((n_init + 1):(n_init + step),:);
    plot(added(:,1),added(:,2),'o', ...
      'Color',[0.000 0.447 0.741], ...
      'MarkerFaceColor',[0.000 0.447 0.741], ...
      'MarkerSize',5,'LineStyle','none');
  end
  if step == 0
    title('initial design');
  else
    title(sprintf('QSI-EEM after %d evaluations',step));
  end
  xlabel('X');
  if k == 1
    ylabel('S');
  end
  xlim([min(x) max(x)]);
  ylim([min(s_regular) max(s_regular)]);
  hold off;

  subplot(2,n_steps,n_steps + k);
  mean_uncertainty = mean(class_uncertainty);
  high_uncertainty = class_uncertainty(:) >= 0.2;
  starts = find(diff([false; high_uncertainty]) == 1);
  stops = find(diff([high_uncertainty; false]) == -1);
  hold on;
  for j = 1:numel(starts)
    x0 = x(starts(j));
    x1 = x(stops(j));
    patch([x0 x1 x1 x0],[-0.02 -0.02 1.02 1.02], ...
      [0.9 0.9 0.9],'EdgeColor','none');
  end
  h_true = plot(x,double(true_proba(:)),'k','LineWidth',1.5);
  h_estimated = plot(x,double(estimated_proba(:)),'Color',[0.000 0.447 0.741], ...
    'LineWidth',1.5);
  h_alpha = plot(x,prm.alpha * ones(size(x)),'Color',[0.8 0.1 0.1], ...
    'LineWidth',1.2);
  ylim([-0.02 1.02]);
  xlim([min(x) max(x)]);
  grid on;
  title(sprintf('set mismatch %.3f, uncertainty %.3f, gray = uncertain', ...
    set_error, mean_uncertainty));
  xlabel('X');
  if k == 1
    ylabel('P(f(x,S) in C)');
    legend([h_true h_estimated h_alpha], {'true','estimated','alpha'}, ...
      'Location','northwest');
  end
  hold off;
end

fig_dir = fullfile(data_dir,'results','figures');
if ~exist(fig_dir,'dir')
  mkdir(fig_dir);
end
paper_write_figure(fig,fig_dir,'qsi_eem_uncertainty_demo');
close(fig);

% --- Comparison of the designs when several algorithms were run
if exist('algorithms','var') && numel(algorithms) > 1
  n_alg = numel(algorithms);
  n_cols = min(n_alg,3);
  n_rows = ceil(n_alg / n_cols);
  fig2 = figure('visible','off');
  clf(fig2);
  set(fig2,'paperpositionmode','auto');
  set(fig2,'position',[100,100,520*n_cols,420*n_rows]);
  for k = 1:n_alg
    prefix = prefixes{k};
    design_k = readmatrix_(fullfile(data_dir,'results','design', ...
      sprintf('doe_%s_%s_%d.csv',prefix,prm.name,run_id)));
    added_k = design_k((n_init + 1):end,:);

    % final GP model of this algorithm, from its saved parameters
    Model = stk_model();
    for m = 1:prm.M
      param_k = readmatrix_(fullfile(data_dir,'results','param', ...
        sprintf('param_%s_%d_%s_%d.csv',prefix,m,prm.name,run_id)));
      cov_k = readmatrix_(fullfile(data_dir,'results','param', ...
        sprintf('cov_%s_%d_%s_%d.csv',prefix,m,prm.name,run_id)));
      Model(m) = stk_model(prm.list_cov{cov_k(end,1)},prm.dim_x + prm.dim_s);
      Model(m).param = param_k(end,:);
    end
    estimated_set_k = get_expected_quantile_set( ...
      Model,df,pts_x,pts_s,design_k,f(design_k),prm.const,prm.alpha);
    set_error_k = lebesgue_deviation(true_set,estimated_set_k);

    subplot(n_rows,n_cols,k);
    p = pcolor(x,s_regular,z_grid);
    set(p,'EdgeColor','none');
    hold on;
    contour(x,s_regular,c_grid,[1 1],'black','LineWidth',1.5);
    plot(initial_design(:,1),initial_design(:,2),'.', ...
      'Color',[0.25 0.25 0.25],'MarkerSize',12);
    plot(added_k(:,1),added_k(:,2),'o', ...
      'Color',[0.000 0.447 0.741], ...
      'MarkerFaceColor',[0.000 0.447 0.741], ...
      'MarkerSize',5,'LineStyle','none');
    title(sprintf('%s: %d evaluations, set mismatch %.3f', ...
      algorithms{k},size(added_k,1),set_error_k));
    xlabel('X');
    if mod(k - 1,n_cols) == 0
      ylabel('S');
    end
    xlim([min(x) max(x)]);
    ylim([min(s_regular) max(s_regular)]);
    hold off;
  end
  paper_write_figure(fig2,fig_dir,'qsi_eem_design_comparison');
  close(fig2);
end

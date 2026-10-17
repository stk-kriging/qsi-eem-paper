% Designs produced by QSI-EEM and QSI-SUR on branin_mod, from the same
% initial design.  Runs both algorithms for config.T steps and writes one
% EPS per method, showing the function, the boundary of f^{-1}(C) (black),
% its posterior-mean estimate (red), and the evaluation points.

close all;

script_dir = fileparts (mfilename ('fullpath'));

data_dir    = fullfile (script_dir, 'outputs');
results_dir = fullfile (data_dir, 'results');
fig_dir     = fullfile (results_dir, 'figures');
paper_prepare_data_dir (data_dir);

case_name = 'branin_mod';
it = 1;
T = 30;         % step shown, = config.T
PTS_DIM = 250;

opts = struct ( ...
    'fig_width_cm',  8.0, ...   % included at 70mm in the paper
    'fig_height_cm', 6.4, ...
    'font_size',     9);

name_list = ["QSI_EEM", "QSI_m"];
name_graphs = ["QSI-EEM", "QSI-SUR"];
fig_names = ["doe_comp_eem", "doe_comp_sur"];

[problem, config] = paper_get_test_case (case_name);
[prm, f, s_trnsf] = problem ();

paper_generate_initial_designs (problem, config, it, data_dir);

paper_set_seed (it);
QSI_EEM (problem, config, it, data_dir);
paper_set_seed (it);
QSI_SUR (problem, config, it, data_dir);

% regular s grid for display, transformed one for the probabilities
xf = stk_sampling_regulargrid (PTS_DIM, prm.dim_x, prm.BOXx);
sf = stk_sampling_regulargrid (PTS_DIM, prm.dim_s, prm.BOXs);
df_disp = double (adapt_set (xf, sf));
zf_disp = f (df_disp);

f_map = reshape (zf_disp, PTS_DIM, PTS_DIM)';
bnd_true = reshape (zf_disp <= prm.const(2,1), PTS_DIM, PTS_DIM)';

df = double (adapt_set (xf, s_trnsf (sf)));
zf = f (df);
true_set = get_true_quantile_set (zf, PTS_DIM, PTS_DIM, prm.alpha, prm.const);

s_lo = min (double (sf));
s_hi = max (double (sf));

for m = 1:numel (name_list)

    name = name_list(m);

    file_doe = fullfile (results_dir, ...
        sprintf ('design/doe_%s_%s_%d.csv', name, prm.name, it));
    file_para = fullfile (results_dir, ...
        sprintf ('param/param_%s_1_%s_%d.csv', name, prm.name, it));
    file_cov = fullfile (results_dir, ...
        sprintf ('param/cov_%s_1_%s_%d.csv', name, prm.name, it));

    doe = readmatrix_ (file_doe);
    para = readmatrix_ (file_para);
    cov_ind = readmatrix_ (file_cov);

    dn = doe(1:config.pts_init + T, :);
    zn = f (dn);

    Model = stk_model (prm.list_cov{cov_ind(T + 1)}, prm.dim_x + prm.dim_s);
    Model.param = para(T + 1, :);

    pred_disp = stk_predict (Model, dn, zn, df_disp);
    bnd_plug = reshape (pred_disp.mean <= prm.const(2,1), PTS_DIM, PTS_DIM)';
    pred = stk_predict (Model, dn, zn, df);
    in_C = (pred.mean >= prm.const(1,1)) & (pred.mean <= prm.const(2,1));
    plug_set = (mean (reshape (in_C, PTS_DIM, PTS_DIM), 2) <= prm.alpha)';

    fig = paper_new_figure (opts);
    ax = paper_new_axes (fig, opts);
    set (ax, 'XGrid', 'off', 'YGrid', 'off', 'Layer', 'top');

    % a raster image, not pcolor: patch seams show through ps2pdf
    image (ax, 'XData', double (xf), 'YData', double (sf), ...
        'CData', f_map, 'CDataMapping', 'scaled');
    colorbar (ax, 'FontSize', opts.font_size);

    contour (ax, double (xf), double (sf)', bnd_true, [1 1], ...
        'LineWidth', 1, 'LineColor', 'black');
    contour (ax, double (xf), double (sf)', bnd_plug, [1 1], ...
        'LineWidth', 1, 'LineColor', 'red');

    plot_set_band (ax, xf, true_set, s_lo + 0.02 * (s_hi - s_lo), 'black');
    plot_set_band (ax, xf, plug_set, s_lo + 0.055 * (s_hi - s_lo), 'red');

    scatter (ax, dn(1:config.pts_init, 1), dn(1:config.pts_init, 2), ...
        12, 'black', 'filled');
    scatter (ax, dn(config.pts_init + 1:end, 1), ...
        dn(config.pts_init + 1:end, 2), 16, 'red', 'filled', ...
        'MarkerEdgeColor', 'black');

    xlim (ax, double (prm.BOXx));
    ylim (ax, double (prm.BOXs));
    xlabel (ax, 'x', 'FontSize', opts.font_size, 'Color', 'black');
    ylabel (ax, 's', 'FontSize', opts.font_size, 'Color', 'black');
    % n includes the initial design, as in the paper.
    title (ax, sprintf ('%s, n = %d', name_graphs(m), ...
        size (dn, 1)), ...
        'FontSize', opts.font_size, 'FontWeight', 'normal');

    paper_save_figure (fig, fig_dir, char (fig_names(m)));
    close (fig);

end

fprintf ('Done.\n');


function plot_set_band (ax, xf, set_ind, y, color)

sel = nan (numel (set_ind), 1);
sel(set_ind == 1) = y;
line (ax, double (xf(:)), sel, 'LineWidth', 3, 'Color', color);

end % function

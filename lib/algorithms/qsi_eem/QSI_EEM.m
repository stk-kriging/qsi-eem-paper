% QSI_EEM sequential design by the EEM criterion, without SMC levels
%
% SYNTAX: QSI_EEM (FUNCT_STRUCT, CONFIG, IT, DATA_DIR, DEMO)
%
%   CONFIG.T points selected on the fixed constraint set PRM.const.
%   Set DEMO to 1 to plot the progression (1D+1D problems only).

% Copyright Notice
%
%    Copyright (C) 2025, 2026 CentraleSupelec
%
%    Author(s): Romain Ait Abdelmalek-Lomenech <romain.ait@centralesupelec.fr>
%               Julien Bect <julien.bect@centralesupelec.fr>
%               Emmanuel Vazquez <emmanuel.vazquez@centralesupelec.fr>

% Copying Permission Statement
%
%    This file is part of qsi-eem-paper (https://github.com/stk-kriging/qsi-eem-paper)
%
%    qsi-eem-paper is free software: you can redistribute it and/or modify it under
%    the terms of the GNU General Public License as published by the Free
%    Software Foundation,  either version 3  of the License, or  (at your
%    option) any later version.
%
%    contrib-qsi is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public
%    License for more details.
%
%    You should  have received a copy  of the GNU  General Public License
%    along with qsi-eem-paper.  If not, see <http://www.gnu.org/licenses/>.

function QSI_EEM (funct_struct, config, it, data_dir, DEMO)

if nargin < 4
    here = fileparts (mfilename ('fullpath'));
    data_dir = fullfile (here, 'data');
end

if nargin < 5
    DEMO = 0;
end

fprintf('Run %d\n', it);

[prm, f, s_trnsf] = funct_struct(); %loading function and parameters
config.batch_size = 1;

if DEMO == 1 && ((prm.dim_x ~= 1) || (prm.dim_s ~= 1))
    error("Invalid problem dimension for demo")
end

if (DEMO ~= 0) && (DEMO ~=1)
    error("Invalid value for variable DEMO.")
end

dim_tot = prm.dim_x+prm.dim_s;

%Initial design
file_grid = sprintf ('doe_init_%s_%d_init.csv', prm.name, it);
di = readmatrix_(fullfile(data_dir, 'doe_init', file_grid));
zi = f(di);

% Create dataframes
dn = stk_dataframe(di);
zn = stk_dataframe(zi);

%stocking parameters
save_param = zeros(config.T+1,dim_tot+1,prm.M);
save_cov = zeros(config.T+1, 1, prm.M);

time = [];

% Estimate and save parameters
Model = stk_model ();
for m = 1:prm.M
    [Model(m), ind_cov] = estim_matern ...
        (dn, zn(:,m), prm.list_cov, config.lognugget);
    save_cov(1,:,m) = ind_cov;
    save_param(1,:,m) = Model(m).param;
end

for t = 1:config.T %loop on steps
    if DEMO == 1
        fprintf('Step %d: evaluating QSI-EEM criterion...\n', t);
    end

    tic

    %sampling points in X x S
    xt = stk_sampling_randomlhs(config.pts_x,prm.dim_x,prm.BOXx);
    st = stk_sampling_randomlhs(config.pts_s,prm.dim_s,prm.BOXs);
    st = s_trnsf(st);
    dt = adapt_set(xt,st);

    %evaluate misclassification probability of \tau(x)
    proba_x = proba_tau(Model, dn , zn, xt, st, prm, config);
    misclass_x = min(proba_x, 1-proba_x);

    [~, sort_misclass_x] = sort(misclass_x,'descend');
    ranking_x = [sort_misclass_x(1)]; %keeping most uncertain point

    %importance sampling on X
    if sum((misclass_x > 0))<=config.keep_x
        ranking_x = sort_misclass_x(1:config.keep_x);
    else
        candidate_set = setdiff(find(misclass_x > 0),ranking_x);
    end

    while size(ranking_x, 1) < config.keep_x

        candidate = randsample(candidate_set, 1, true, misclass_x(candidate_set));

        candidate_set = setdiff(candidate_set, candidate); %enforce without replacement constraint
        ranking_x = [ranking_x; candidate];
    end

    %dt = adapt_set(xt(ranking_x,:), st);
    xt = xt(ranking_x,:);
    dt = adapt_set(xt, st);

    % Define IS weights
    IS = 1 ./ misclass_x(ranking_x);
    IS(~isfinite(IS)) = 0;
    if sum(IS) > 0
        IS = IS/sum(IS);
    else
        IS = ones(size(IS))/numel(IS);
    end

    % plug-in estimator on the current grid
    mu_n = predict_mean_var(Model, dn, zn, dt, prm);
    estim_current = get_estim(mu_n, prm.const, prm, config);

    newpt = get_point_EEM(Model, estim_current, dn, zn, xt, st, prm.const, IS, prm, config, 0);

    if DEMO == 1
        fprintf('Evaluating function at (%g , %g)\n', newpt(1), newpt(2));
    end

    time = [time, toc];

    dn = stk_dataframe([dn;newpt]);
    zn = stk_dataframe([zn;f(newpt)]);

    if DEMO == 1
        disp('Updating GP parameters...')
    end

    if mod(t, config.estim_param_steps) == 0
        for m = 1:prm.M
            [Model(m), ind_cov] = estim_matern ...
                (dn, zn(:,m), prm.list_cov, config.lognugget);
            save_cov(t+1,:,m) = ind_cov;
            save_param(t+1,:,m) = Model(m).param;
        end
    else
        save_cov(t+1,:,:) = save_cov(t,:,:);
        save_param(t+1,:,:) = save_param(t,:,:);
    end

    % --- Save results for this step ---
    write_results(data_dir, prm, it, dn, save_param, save_cov, time);

    if DEMO == 1

        PLOT_PERIOD = 3;

        if (mod (t, PLOT_PERIOD) == 0) && (t >= PLOT_PERIOD)
            fprintf("Generating graphs (%d points)...\n", t);

            if t == PLOT_PERIOD

                h0 = get (0);
                w = h0.MonitorPositions(3) / 12;
                h = h0.MonitorPositions(4);
                y = 0.8 * h;

                figs(1) = figure ('Position', [1.50*w y   w 0.4*h]);
                figs(2) = figure ('Position', [2.75*w y 4*w 0.4*h]);
                figs(3) = figure ('Position', [7.00*w y 4*w 0.4*h]);

                A = stk_sampling_regulargrid(1000, 1, prm.BOXs);
                A = sort(double(A));
                [~, pdf] = s_trnsf(A);

                figure (figs(1));
                hold on
                area(A, pdf, 'FaceAlpha',0.2);
                ylabel("\bfDens. func.")
                xticks([])
                camroll(90);
                drawnow ();
            end

            make_graphs_ (figs, data_dir, funct_struct, config, ...
                sprintf("QSI_%s", config.critName), sprintf("%d steps", t), it, t, 0);
            disp("PAUSED: press any key to continue.")
            pause()
        end

    end

end

end


function write_results(data_dir, prm, it, dn, save_param, save_cov, time)
% Write the current design, GP parameters and timings to disk.

    filename = sprintf ('doe_QSI_EEM_%s_%d.csv', prm.name, it);
    writematrix_ (double (dn), fullfile (data_dir, 'results/design', filename));

    for m = 1:prm.M
        filename = sprintf ('param_QSI_EEM_%d_%s_%d.csv', m, prm.name, it);
        writematrix_ (save_param(:,:,m), fullfile (data_dir, 'results/param', filename));

        filename = sprintf ('cov_QSI_EEM_%d_%s_%d.csv', m, prm.name, it);
        writematrix_ (save_cov(:,:,m), fullfile (data_dir, 'results/param', filename));
    end

    filename = sprintf ('time_QSI_EEM_%s_%d.csv', prm.name, it);
    writematrix_ (time, fullfile (data_dir, 'results/time', filename));

end % function


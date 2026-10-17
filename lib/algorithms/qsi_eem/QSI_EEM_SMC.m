% QSI_EEM_SMC sequential design for small quantile sets (QSI-EEM + SMC)
%
% SYNTAX: QSI_EEM_SMC (FUNCT_STRUCT, CONFIG, IT, BATCH_SIZE, DATA_DIR, MUTE)
%
%   Tracks the target set through a sequence of relaxed constraint sets,
%   selecting batches of evaluation points by the EEM criterion.  IT names
%   the initial design and suffixes the outputs, written under
%   DATA_DIR/results after each step.  Set MUTE to 1 for quiet runs.

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

function QSI_EEM_SMC (funct_struct, config, it, batch_size, data_dir, mute)

if nargin < 5
    here = fileparts (mfilename ('fullpath'));
    data_dir = fullfile (here, 'data');
end

if nargin < 6
    mute = 0;
end

fprintf('Run %d\n', it);

[prm, f, s_trnsf] = funct_struct(); %loading function and parameters
config.batch_size = batch_size;
config.keep_x = config.pts_x; % no subsampling, whatever pts_x is

dim_tot = prm.dim_x+prm.dim_s;

%Initial design
file_grid = sprintf ('doe_init_%s_%d_init.csv', prm.name, it);
di = readmatrix_(fullfile(data_dir, 'doe_init', file_grid));
zi = f(di);

% Create dataframes
dn = stk_dataframe(di);
zn = stk_dataframe(zi);

%stocking parameters
save_param = zeros(1,dim_tot+1,prm.M);
save_cov = zeros(1, 1, prm.M);

time = [];

% Estimate and save parameters
Model = stk_model ();
for m = 1:prm.M
    [Model(m), ind_cov] = estim_matern ...
        (dn, zn(:,m), prm.list_cov, config.lognugget);
    save_cov(1,:,m) = ind_cov;
    save_param(1,:,m) = Model(m).param;
end

st = stk_sampling_sobol(config.pts_s,prm.dim_s,prm.BOXs);
st = s_trnsf(st);

smc = struct();

final_stage = 0;
ending = 0;
resample = 0;
t = 0;

while ending ~= 1
    t = t+1;
    if ~mute
        fprintf('Step %d\n', t);
    end
    tic

    if t == 1
        smc = smc_init(Model, dn, zn, st, prm, config);
    elseif resample == 1
        smc = smc_next_level(smc, Model, dn, zn, st, prm, config, t);
        resample = 0;
    end

    if all(smc.const(2,:) == prm.const(2,:))
        final_stage = 1;
    end

    if ~mute
        if final_stage
            level_str = sprintf('%d (final)', numel(smc.const_list));
        else
            level_str = sprintf('%d', numel(smc.const_list));
        end
        fprintf('  Level %s: constraint set %s, beta = %.3g, sigma = %s\n', ...
            level_str, mat2str(double(smc.const), 5), smc.beta, mat2str(smc.sigma, 3));
    end

    % --- Evaluate the EEM criterion and pick the next point ---
    if ~mute
        fprintf('  Evaluating EEM criterion\n');
    end
    % uniform IS weights
    IS = 1/size(smc.xt,1)*ones(size(smc.xt,1),1);

    newpt = get_point_EEM(Model, smc.estim_current, dn, zn, smc.xt, st, smc.const, IS, prm, config, 1);

    time = [time, toc];

    % --- Evaluate f at the new point and re-estimate the GP models ---
    if ~mute
        fprintf('  EEM criterion time: %.1f s\n', time(end));
        n = size(dn, 1);
        for i = 1:size(newpt, 1)
            fprintf('  Evaluating f at point %d: %s\n', ...
                n+i, mat2str(double(newpt(i,:)), 5));
        end
    end
    dn = stk_dataframe([dn;newpt]);
    zn = stk_dataframe([zn;f(newpt)]);

    if mod(t, config.estim_param_steps) == 0
        if ~mute
            fprintf('  Estimating GP parameters\n');
        end
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

    [mu_n, var_n] = predict_mean_var(Model, dn, zn, smc.dt, prm);
    smc.mu_n = mu_n;  smc.var_n = var_n;  % reused at the next level change
    smc.estim_current = get_estim(mu_n, smc.const, prm, config);
    estim_safe = get_estim_safe(mu_n, var_n, smc.beta, smc.const, prm.alpha, config.pts_x, config.pts_s);

    if ~mute
        fprintf('  Particles in estimated set: %.3f (safe estimator: %.3f)\n', ...
            mean(smc.estim_current), mean(estim_safe));
    end

    % --- Stopping criterion: tighten the constraint, or stop at final stage ---
    if stop_criterion(final_stage, Model, dn, zn, smc.dt, smc.estim_current, mu_n, IS, smc.const, prm, config)
        if ~final_stage
            resample = 1;
            if ~mute
                fprintf('  Stopping criterion met: tightening the constraint set\n');
            end
        else
            ending = 1;
            if ~mute
                fprintf('  Stopping criterion met at final stage: ending\n');
            end
        end
    end

    % --- Restart criterion: roll back to an earlier constraint set ---
    if (mean(estim_safe) > 2*config.rho) || (mean(estim_safe) < config.rho/2)
        if ~mute
            fprintf('  Restart criterion triggered\n');
        end
        [smc, restarted] = smc_restart(smc, Model, dn, zn, st, IS, prm, config, t, mute, estim_safe);
        if restarted
            final_stage = 0;
            resample = 0;
            ending = 0;
        else
            smc.level_log(end+1,:) = [t, 5, numel(smc.const_list), ...
                smc.thresh, smc.beta, mean(smc.estim_current), mean(estim_safe)];
        end
    end

    % --- Save results for this step ---
    write_smc_results(data_dir, prm, it, batch_size, dn, ...
        save_param, save_cov, time, smc.level_log);

    if isfield(config, 'max_steps') && (t >= config.max_steps)
        if ~mute && ~ending
            fprintf('  Maximum number of steps reached\n');
        end
        ending = 1;
    end

end

end


function smc = smc_init(Model, dn, zn, st, prm, config)
% initial population and constraint set

    smc.xt = stk_sampling_maximinlhs(config.pts_x, prm.dim_x, prm.BOXx);
    smc.dt = adapt_set(smc.xt, st);

    [smc.const, smc.thresh, mu_n] = init_set(Model, smc.dt, dn, zn, prm, config);
    smc.sigma = config.sigma_init;
    smc.beta = 1/2;
    smc.estim_current = get_estim(mu_n, smc.const, prm, config);

    smc.x_list = {smc.xt};
    smc.beta_list = smc.beta;
    smc.const_list = smc.thresh;

    % level_log 
    % rows: [step, event, level, thresh, beta, prop, prop_safe]
    % events: 0 init, 1 new level, 2 restart-accept, 3 restart-rebuild, 4 restart-rebuild-init,
    % 5 restart-dismissed, 6 restart-refresh
    smc.level_log = [1, 0, 1, smc.thresh, smc.beta, mean(smc.estim_current), NaN];

end


function smc = smc_next_level(smc, Model, dn, zn, st, prm, config, t)
% resample and move the particles, then tighten the constraint

    smc.beta = calibrate_kappa(Model, smc.dt, dn, zn, smc.const, config, prm, ...
        smc.mu_n, smc.var_n);
    smc.beta_list = [smc.beta_list, smc.beta];

    estim_safe = get_estim_safe(smc.mu_n, smc.var_n, smc.beta, smc.const, prm.alpha, config.pts_x, config.pts_s);
    smc.xt = resample_xt(smc.xt, estim_safe, config);

    [smc.xt, smc.sigma] = MH_walk(Model, smc.xt, st, dn, zn, smc.const, smc.sigma, smc.beta, prm, config);
    smc.x_list{end+1} = smc.xt;
    smc.dt = adapt_set(smc.xt, st);

    [smc.const, smc.thresh, mu_n] = new_const(Model, smc.dt, dn, zn, smc.thresh, smc.beta, prm, config);
    smc.const_list = [smc.const_list; smc.thresh];
    smc.estim_current = get_estim(mu_n, smc.const, prm, config);

    smc.level_log(end+1,:) = [t, 1, numel(smc.const_list), smc.thresh, ...
        smc.beta, mean(smc.estim_current), NaN];

end


function [smc, restarted] = smc_restart(smc, Model, dn, zn, st, IS, prm, config, t, mute, estim_safe_in)
% roll back to the first level that is no longer valid; rebuild the
% last one if none is found (or keep it if it is the final set,
% refreshing its particles when too few remain inside)

    restarted = true;
    smc_in = smc;
    nb_levels = length(smc.const_list);

    for u = 1:nb_levels
        smc.xt = smc.x_list{u};
        smc.dt = adapt_set(smc.xt, st);
        smc.const = threshold_to_const(prm, smc.const_list(u));

        if u == nb_levels
            if all(smc.const(2,:) == prm.const(2,:))
                if mean(estim_safe_in) < config.rho/2
                    % population lost the final set: resample what remains and move
                    if ~mute
                        fprintf('  Restart: refreshing the particles of the final constraint set\n');
                    end
                    smc = smc_in;
                    if any(estim_safe_in)
                        smc.xt = resample_xt(smc.xt, estim_safe_in, config);
                    end
                    [smc.xt, smc.sigma] = MH_walk(Model, smc.xt, st, dn, zn, smc.const, smc.sigma, smc.beta, prm, config);
                    smc.x_list{end} = smc.xt;
                    smc.dt = adapt_set(smc.xt, st);
                    mu_n = predict_mean_var(Model, dn, zn, smc.dt, prm);
                    smc.estim_current = get_estim(mu_n, smc.const, prm, config);
                    smc.level_log(end+1,:) = [t, 6, numel(smc.const_list), ...
                        smc.thresh, smc.beta, mean(smc.estim_current), mean(estim_safe_in)];
                    return
                end
                if ~mute
                    fprintf('  Restart: no problematic level, keeping the final constraint set\n');
                end
                smc = smc_in;
                restarted = false;
                return
            end
            if ~mute
                if u == 1
                    fprintf('  Restart: rebuilding the initial constraint set\n');
                else
                    fprintf('  Restart: rebuilding level %d from level %d\n', u, u-1);
                end
            end
            if u == 1
                % no earlier level: rebuild the initial set
                event = 4;
                [smc.const, smc.thresh, mu_n] = init_set(Model, smc.dt, dn, zn, prm, config);
                smc.beta = 1/2;
                smc.beta_list = smc.beta;
                smc.const_list = smc.thresh;
            else
                event = 3;
                smc.const_list = smc.const_list(1:u-1);
                smc.beta_list = smc.beta_list(1:u-1);
                smc.const = threshold_to_const(prm, smc.const_list(u-1));
                smc.thresh = smc.const_list(u-1);

                smc.beta = calibrate_kappa(Model, smc.dt, dn, zn, smc.const, config, prm);
                smc.beta_list = [smc.beta_list, smc.beta];

                [smc.const, smc.thresh, mu_n] = new_const(Model, smc.dt, dn, zn, smc.thresh, smc.beta, prm, config);
                smc.const_list = [smc.const_list; smc.thresh];
            end
            smc.estim_current = get_estim(mu_n, smc.const, prm, config);
            smc.level_log(end+1,:) = [t, event, numel(smc.const_list), ...
                smc.thresh, smc.beta, mean(smc.estim_current), NaN];

        else
            [mu_n, var_n] = predict_mean_var(Model, dn, zn, smc.dt, prm);
            smc.beta = calibrate_kappa(Model, smc.dt, dn, zn, smc.const, config, prm, mu_n, var_n);
            smc.estim_current = get_estim(mu_n, smc.const, prm, config);
            estim_safe = get_estim_safe(mu_n, var_n, smc.beta, smc.const, prm.alpha, config.pts_x, config.pts_s);

            if (mean(estim_safe) > 2*config.rho) || (mean(estim_safe) < config.rho/2) || ~stop_criterion(0, Model, dn, zn, smc.dt, smc.estim_current, mu_n, IS, smc.const, prm, config)
                if ~mute
                    fprintf('  Restart: back to level %d\n', u);
                end
                smc.const_list = smc.const_list(1:u);
                smc.x_list = smc.x_list(1:u);
                smc.beta_list = [smc.beta_list(1:u-1), smc.beta];
                smc.thresh = smc.const_list(u);
                smc.level_log(end+1,:) = [t, 2, u, smc.thresh, smc.beta, ...
                    mean(smc.estim_current), mean(estim_safe)];
                break
            end
        end
    end

end


function write_smc_results(data_dir, prm, it, batch_size, dn, ...
        save_param, save_cov, time, level_log)
% Write the current design, GP parameters, timings and SMC history to disk.

    filename = sprintf ('doe_QSI_EEM_smc_%s_%d_%d.csv', prm.name, it, batch_size);
    writematrix_ (double (dn), fullfile (data_dir, 'results/design', filename));

    for m = 1:prm.M
        filename = sprintf ('param_QSI_EEM_smc_%d_%s_%d_%d.csv', m, prm.name, it, batch_size);
        writematrix_ (save_param(:,:,m), fullfile (data_dir, 'results/param', filename));

        filename = sprintf ('cov_QSI_EEM_smc_%d_%s_%d_%d.csv', m, prm.name, it, batch_size);
        writematrix_ (save_cov(:,:,m), fullfile (data_dir, 'results/param', filename));
    end

    filename = sprintf ('time_QSI_EEM_smc_%s_%d_%d.csv', prm.name, it, batch_size);
    writematrix_ (time, fullfile (data_dir, 'results/time', filename));

    filename = sprintf ('log_QSI_EEM_smc_%s_%d_%d.csv', prm.name, it, batch_size);
    writematrix_ (double (level_log), fullfile (data_dir, 'results/smc_logs', filename));

end % function


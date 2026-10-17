% PAPER_EXTRACT_PRECISION_SMALL_VARIABLE_PARAMS is extract_precision_small for result files suffixed
% with __<parameter>_<value> by paper_rename_files_with_param. VALUE is that
% suffix, as a string.

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
%    qsi-eem-paper is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public
%    License for more details.
%
%    You should  have received a copy  of the GNU  General Public License
%    along with qsi-eem-paper.  If not, see <http://www.gnu.org/licenses/>.
function [deviation, nb_batches] = paper_extract_precision_small_variable_params ...
    (funct_struct, config, algo_name, parameter, value, it, batch_size, data_dir)

fprintf('Run %d\n', it);

[prm, f, s_trnsf] = funct_struct();

st = stk_sampling_sobol(config.pts_s_deter, prm.dim_s, prm.BOXs);
st = s_trnsf(st);

% Metamodel

doe_file = sprintf("doe_%s_%s_%d_%d__%s_%s.csv", ...
    algo_name, prm.name, it, batch_size, parameter, value);
dn = readmatrix_(fullfile(data_dir, "results/design", doe_file));
zn = f(dn);

Model = stk_model();
for m = 1:prm.M
    Model(m) = estim_matern(dn, zn(:,m), prm.list_cov, config.lognugget);
end

f_pred = @(d) posterior_means(Model, dn, zn, d);

% Estimation of vol(Gamma_pred \Delta Gamma) / vol(Gamma)

[x_true, vol_levels_true] = splitting_pass(f, st, prm, config, 0);

vol_true_set = vol_levels_true(end);

if isempty(x_true) || vol_true_set == 0
    warning('Run %d: splitting failed to reach the true quantile set.', it);
    deviation = NaN;
    nb_batches = (size(dn, 1) - config.pts_init) / batch_size;
    return
end

[x_pred, vol_levels_pred] = ...
    splitting_pass(f_pred, st, prm, config, eps(vol_true_set));

dt = adapt_set(x_true, st);
in_true = get_estim_deter(eval_product(f, x_true, st), prm.const, prm, config);
in_pred = get_estim_deter(f_pred(dt), prm.const, prm, config);
vol_missed = outer_volume(vol_levels_true) * mean(in_true & ~in_pred);

if isempty(x_pred)
    vol_spurious = 0;
else
    dt = adapt_set(x_pred, st);
    in_true = get_estim_deter(eval_product(f, x_pred, st), prm.const, prm, config);
    in_pred = get_estim_deter(f_pred(dt), prm.const, prm, config);
    vol_spurious = outer_volume(vol_levels_pred) * mean(in_pred & ~in_true);
end

deviation = (vol_missed + vol_spurious) / vol_true_set;
nb_batches = (size(dn, 1) - config.pts_init) / batch_size;

end % function

function [xt, vol_levels] = splitting_pass(g, st, prm, config, volume_floor)

xt = stk_sampling_maximinlhs(config.pts_x_deter, prm.dim_x, prm.BOXx);

vol_levels = [];
sigma = config.sigma_init;
thresh = NaN;
const = prm.const;
in_current_level = [];
final_level_reached = false;

while ~final_level_reached

    if isempty(vol_levels)
        zt = eval_product(g, xt, st);
        [const, thresh] = deterministic_set_const_init(zt, prm, config);
    else
        xt = resample_xt_deter(xt, in_current_level, config);
        [xt, sigma] = MH_walk_deter(g, xt, st, const, sigma, prm, config);
        zt = eval_product(g, xt, st);
        [const, thresh] = deterministic_set_const(zt, thresh, prm, config);
    end

    final_level_reached = (thresh == final_threshold (prm));

    in_current_level = get_estim_deter(zt, const, prm, config);

    if isempty(vol_levels)
        vol_levels = mean(in_current_level);
    else
        vol_levels = [vol_levels; vol_levels(end) * mean(in_current_level)];
    end

    if (~any(in_current_level)) && (~final_level_reached)
        xt = [];
        return
    end

    if (~final_level_reached) && (vol_levels(end) <= volume_floor)
        xt = [];
        return
    end

end

end % function

function vol = outer_volume(vol_levels)

if length(vol_levels) > 1
    vol = vol_levels(end-1);
else
    vol = 1;
end

end % function

function z = posterior_means(Model, dn, zn, d)

% Posterior means through precomputed kriging weights (constant trend).

block_size = 1e5;
n = size(d, 1);
z = zeros(n, size(zn, 2));

xi = double(dn);
for m = 1:size(zn, 2)
    [K, P] = stk_make_matcov(Model(m), xi);
    if size(P, 2) ~= 1 || any(P ~= P(1))
        error('posterior_means expects a constant trend.');
    end
    zi = double(zn(:,m));
    L = chol((K + K')/2, 'lower');
    beta = (P' * (L' \ (L \ P))) \ (P' * (L' \ (L \ zi)));
    w{m} = L' \ (L \ (zi - P * beta));        %#ok<AGROW>
    trend(m) = P(1) * beta;                   %#ok<AGROW>
end

for i0 = 1:block_size:n
    i1 = min(i0 + block_size - 1, n);
    for m = 1:size(zn, 2)
        Kt = stk_make_matcov(Model(m), d(i0:i1, :), xi);
        z(i0:i1, m) = trend(m) + Kt * w{m};
    end
end

end % function

% DETERMINISTIC_SUBSET estimates the volume of the quantile set
%
% SYNTAX: [VOL, THRESH_LIST, PROP_LIST] = deterministic_subset ...
%             (CASE_NAME, RUN_ID, CONFIG_OVERRIDES)
%
%   estimates lambda(Gamma) / lambda(X) by subset simulation, for the test
%   case CASE_NAME (as understood by PAPER_GET_TEST_CASE).  RUN_ID seeds the
%   run and defaults to 1.  CONFIG_OVERRIDES is an optional structure whose
%   fields overwrite those of the configuration of the test case.
%
%   VOL is the estimated volume, THRESH_LIST and PROP_LIST the threshold and
%   the proportion of particles retained at each level.
%
%   Example:
%
%      vol = deterministic_subset ('rotor37');
%
% See also: paper_get_test_case, extract_precision_small

% Copyright Notice
%
%    Copyright (C) 2026 CentraleSupelec
%
%    Author(s): Romain Ait Abdelmalek-Lomenech <romain.ait@centralesupelec.fr>

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

function [vol, thresh_list, prop_list] = deterministic_subset ...
    (case_name, run_id, config_overrides)

if nargin < 2
    run_id = 1;
end

[problem, config] = paper_get_test_case (case_name);
config.pts_x_deter = 10000;

if nargin > 2
    names = fieldnames (config_overrides);
    for i = 1:(numel (names))
        config.(names{i}) = config_overrides.(names{i});
    end
end

[prm, f, s_trnsf] = problem();

paper_set_seed(run_id);

st = stk_sampling_sobol(config.pts_s_deter, prm.dim_s, prm.BOXs);
st = s_trnsf(st);

xt = stk_sampling_maximinlhs(config.pts_x_deter, prm.dim_x, prm.BOXx);

thresh_list = [];
prop_list = [];
sigma = config.sigma_init;
thresh = NaN;
const = prm.const;
in_current_level = [];
final_level_reached = false;

while ~final_level_reached

    if isempty(prop_list)
        zt = eval_product(f, xt, st);
        [const, thresh] = deterministic_set_const_init(zt, prm, config);
    else
        xt = resample_xt_deter(xt, in_current_level, config);
        [xt, sigma] = MH_walk_deter(f, xt, st, const, sigma, prm, config);
        zt = eval_product(f, xt, st);
        [const, thresh] = deterministic_set_const(zt, thresh, prm, config);
    end

    final_level_reached = (thresh == final_threshold (prm));

    in_current_level = get_estim_deter(zt, const, prm, config);
    thresh_list = [thresh_list; thresh];
    prop_list = [prop_list; mean(in_current_level)];

    fprintf('Level %d: threshold %.6g, proportion %.3f, volume %.3g\n', ...
        numel(prop_list), thresh, prop_list(end), prod(prop_list));

    if (~any(in_current_level)) && (~final_level_reached)
        error('splitting stalled before reaching the quantile set');
    end

end

vol = prod(prop_list);
fprintf('Estimated volume of the quantile set: %.4g (%d levels)\n', ...
    vol, numel(thresh_list));

end % function

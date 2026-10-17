% CALIBRATE_KAPPA calibrates the relaxation parameter beta

% Copyright Notice
%
%    Copyright (C) 2025 CentraleSupelec
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

function beta_act = calibrate_kappa (Model, dt, dn, zn, const, config, prm, mu_n, var_n)

    % no relaxation requested
    if config.kappa <= 0
        beta_act = 1/2;
        return
    end

    is_valid = 0;

    beta_0 = 1/2;
    beta_1 = 1;
    beta_act = 0.75;

    if nargin < 9
        [mu_n, var_n] = predict_mean_var(Model, dn, zn, dt, prm);
    end

    estim = get_estim(mu_n, const, prm, config);

    % empty estimated set: keep the default beta
    if mean(estim) == 0
        beta_act = 1/2;
        return
    end

    iterations = 0;
    max_iterations = 60;

    while is_valid ~= 1
    iterations = iterations + 1;
    estim_safe = get_estim_safe(mu_n, var_n, beta_act, const, prm.alpha, config.pts_x, config.pts_s);
    if (beta_act/(1/2) < 1.0001) || (beta_act > 0.9999)
        beta_act = 1/2;
        is_valid = 1;
    elseif mean(estim_safe)/mean(estim) >= 1 + config.kappa + 0.025
        beta_1 = beta_act;
        beta_act = (beta_0+beta_act)/2;
    elseif mean(estim_safe)/mean(estim) <= 1 + config.kappa - 0.025
        beta_0 = beta_act;
        beta_act = (beta_1+beta_act)/2;

    else
        is_valid = 1;
    end

    % the target ratio may be unreachable
    if iterations >= max_iterations
        is_valid = 1;
    end

    end
end % function


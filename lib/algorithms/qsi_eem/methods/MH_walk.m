% MH_WALK performs a few steps of Metropolis-Hastings on the safe density

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

function [xt, sigma] = MH_walk (Model, ...
    xt, st, dn, zn, const, sigma, beta, prm, config)

    for j = 1:config.MH_steps

        xt_prop = xt + sigma.*randn(size(xt,1),prm.dim_x);

        % only predict for proposals inside the box
        in_box = ismember(xt_prop, prm.BOXx);
        is_valid = false(size(xt,1), 1);

        if any(in_box)
            dt = adapt_set(xt_prop(in_box,:), st);
            [mu_n, var_n] = predict_mean_var(Model, dn, zn, dt, prm);
            is_valid(in_box) = get_estim_safe(mu_n, var_n, beta, const, prm.alpha, sum(in_box), config.pts_s);
        end

        xt = is_valid.*xt_prop + (1-is_valid).*xt;

        if mean(is_valid) > config.MH_accept
            sigma = exp(log(sigma) + config.MH_delta);
        else
            sigma = exp(log(sigma) - config.MH_delta);
        end
    end

end % function


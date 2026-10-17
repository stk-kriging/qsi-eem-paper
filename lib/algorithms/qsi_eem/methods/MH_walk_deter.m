% MH_WALK_DETER Metropolis walk of the particles, on simulator values

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

function [xt, sigma] = MH_walk_deter (f, xt, st, const, sigma, prm, config)

    for j = 1:config.MH_steps
        xt_prop = xt + sigma.*randn(size(xt,1),prm.dim_x);

        % only evaluate the proposals inside the box, as in MH_walk
        in_box = ismember(xt_prop, prm.BOXx);
        estim_safe = false(size(xt,1), 1);
        if any(in_box)
            zt = eval_product(f, xt_prop(in_box,:), st);
            estim_safe(in_box) = get_estim_deter(zt, const, prm, config);
        end
        is_valid = (in_box & estim_safe);

        xt = is_valid.*xt_prop + (1-is_valid).*xt;

        if mean(is_valid) > config.MH_accept
            sigma = exp(log(sigma) + config.MH_delta);
        else
            sigma = exp(log(sigma) - config.MH_delta);
        end

    end

end % function


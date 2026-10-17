% GET_ESTIM membership of each particle in the estimated quantile set

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

function [estim, temp_box] = get_estim (mu_n, const, prm, config)

    eval = zeros(size(mu_n,1),1);
    temp_box = const;

    for m = 1:prm.M  % assume outputs are independent
        mu = mu_n(:,m);

        if nargout > 1
            temp_box(:,m) = const(:,m);
        end

        if const(1,m) > const(2,m)
            eval_temp = (const(1,m) < mu) | (const(2,m) > mu);
        else
            eval_temp = (const(1,m) < mu) & (const(2,m) > mu);
        end

        eval = (eval | eval_temp);
    end

    % one row per particle of the x population
    eval = reshape(eval, size(mu_n,1)/config.pts_s, config.pts_s);

    proba = sum(eval,2) / config.pts_s;
    estim = (proba <= prm.alpha);

end % function


% GET_PROBA_DETER proportion of the S sample in the critical region
%
% PROBA = GET_PROBA_DETER (MU_N, CONST, PRM, CONFIG) returns, for each x
% of a product design with x varying fastest (see ADAPT_SET), the
% proportion of the CONFIG.PTS_S_DETER points of S where the outputs
% MU_N lie in the critical region defined by CONST.  Same rule as
% GET_ESTIM_DETER, which returns the thresholded set PROBA <= PRM.ALPHA.

% Copyright Notice
%
%    Copyright (C) 2026 CentraleSupelec
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

function proba = get_proba_deter (mu_n, const, prm, config)

    eval = false (size (mu_n, 1), 1);

    for m = 1:prm.M  % assume outputs are independent
        mu = mu_n(:,m);

        if const(1,m) > const(2,m)
            eval_temp = (const(1,m) < mu) | (const(2,m) > mu);
        else
            eval_temp = (const(1,m) < mu) & (const(2,m) > mu);
        end

        eval = (eval | eval_temp);
    end

    eval = reshape (eval, size (mu_n, 1) / config.pts_s_deter, config.pts_s_deter);

    proba = sum (eval, 2) / config.pts_s_deter;

end % function

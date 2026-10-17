% PREDICT_MEAN_VAR posterior means (and variances) of the M outputs

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

function [mu_n, var_n] = predict_mean_var (Model, dn, zn, dt, prm)

    mu_n = zeros(size(dt,1), prm.M);
    want_var = (nargout > 1);
    if want_var
        var_n = zeros(size(dt,1), prm.M);
    end

    for m = 1:prm.M
        pred = stk_predict(Model(m), dn, zn(:,m), dt);
        mu_n(:, m) = pred.mean;
        if want_var
            var_n(:,m) = pred.var;
        end
    end

end % function


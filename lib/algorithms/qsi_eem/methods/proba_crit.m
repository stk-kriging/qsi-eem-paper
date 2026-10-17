% PROBA_CRIT probability of falling between the constraint bounds

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

function proba = proba_crit (Model, dn, zn, dt, const, prm)

    proba = 1;

    for m = 1:prm.M
        pred = stk_predict(Model(m),dn,zn(:,m),dt);
        proba = proba .* abs(stk_distrib_normal_cdf(const(2,m),pred.mean,sqrt(pred.var)) - stk_distrib_normal_cdf(const(1,m),pred.mean,sqrt(pred.var)))';
    end

end % function


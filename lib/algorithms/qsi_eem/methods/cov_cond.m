% COV_COND posterior covariance between two sets of points

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

function post_cov = cov_cond (model, dn, zn, dt, dt_cand)

prior_model = model;
post_model = stk_model_gpposterior (model, dn, zn);
kreq1 = stk_make_kreq (post_model, dt);
kreq2 = stk_make_kreq (post_model, dt_cand);

prior_cov = stk_make_matcov (prior_model, dt, dt_cand);

delta_cov = (get (kreq1, 'lambda_mu'))' * (get (kreq2, 'RS'));
post_cov = prior_cov - delta_cov;

end % function


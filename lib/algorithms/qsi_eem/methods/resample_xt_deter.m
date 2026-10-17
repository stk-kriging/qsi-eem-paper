% RESAMPLE_XT_DETER residual resampling of the particles in the current level
%
%   Same as RESAMPLE_XT, with CONFIG.pts_x_deter particles.

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

function xt = resample_xt_deter (xt, in_set, config)

    idx = find (double (in_set));
    idx = idx(:);
    K = numel (idx);
    M = config.pts_x_deter;

    if K == 0
        error ('resample_xt_deter: no particle in the set');
    end

    % deterministic part
    n_det = floor (M / K);
    resampling = repmat (idx, n_det, 1);

    % residual part
    R = M - K * n_det;
    if R > 0
        resampling = [resampling; idx(randi (K, R, 1))];
    end

    xt = xt(resampling, :);

end % function

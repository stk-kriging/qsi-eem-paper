% STOP_CRITERION decides whether the current constraint set is resolved

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

function stop = stop_criterion (final_stage, Model, ...
    dn, zn, dt, estim_current, mu_n, IS, const, prm, config)

    ranking_xs = get_candidates(Model, dn, zn, dt, const, prm, config, 1);
    pt_temp = dt(ranking_xs,:);
    pt = [];

    % keep the candidates not already in the design
    for i = 1:size(pt_temp,1)
        if size(unique(double([dn; pt; pt_temp(i,:)]), "rows"), 1) == size(dn,1)+size(pt,1)+1
            pt = [pt; pt_temp(i,:)];
        end
    end

    if mean(estim_current) == 0
        stop = 0;
        return
    end

    potential = modif_mu(pt, [], Model, dn, zn, dt, estim_current, mu_n, IS, const, prm, config)/mean(estim_current);

    if final_stage && (potential < config.stop_final)
        stop = 1;
    elseif ~final_stage && (potential < config.stop)
        stop = 1;
    else
        stop = 0;
    end

end % function


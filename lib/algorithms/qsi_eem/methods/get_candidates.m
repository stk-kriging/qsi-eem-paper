% GET_CANDIDATES ranks candidate points by misclassification probability

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

function ranking_xs = get_candidates ...
    (Model, dn, zn, dt, const, prm, config, stop)

    proba_xs = proba_crit(Model, dn, zn, dt, const, prm);
    misclass_xs = min(proba_xs,1-proba_xs)';

    % more candidates when checking the stopping criterion
    if stop == 0
        target = config.keep_xs;
    else
        target = config.size_stop_test;
    end

    [~, sort_misclass_xs] = sort(misclass_xs, 'descend');
    positive = find(misclass_xs > 0);

    if numel(positive) <= target
        ranking_xs = sort_misclass_xs(1:target);
    else
        % weighted sampling without replacement, keeping the top point
        top = sort_misclass_xs(1);
        pool = setdiff(positive, top);
        keys = rand(numel(pool), 1).^(1./misclass_xs(pool));
        [~, order] = sort(keys, 'descend');
        ranking_xs = [top; pool(order(1:target-1))];
    end
end % function


% GET_POINT_EEM selects the next batch of evaluation points

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

function pts_to_evaluate = get_point_EEM(Model, estim_current, dn, zn, xt, st, const, IS, prm, config, optimize)

    old_pts = [];

    % precomputations for modif_mu
    precomp = struct();
    for m = 1:prm.M
        precomp.post_model{m} = stk_model_gpposterior(Model(m), dn, zn(:,m));
    end

    for b = 1:config.batch_size
        idx = randsample(1:size(xt,1), config.keep_x, false);
        estim_current_temp = estim_current(idx);
        xt_temp = xt(idx, :);
        dt = adapt_set(xt_temp, st);

        % IS weights of the subsample
        IS_temp = IS(idx);
        if sum(IS_temp) > 0
            IS_temp = IS_temp/sum(IS_temp);
        else
            IS_temp = ones(numel(idx), 1)/numel(idx);
        end

        mu_n = predict_mean_var(Model, dn, zn, dt, prm);

        for m = 1:prm.M
            kreq_dt = stk_make_kreq(precomp.post_model{m}, dt);
            precomp.lambda_mu_dt{m} = (get(kreq_dt, 'lambda_mu'))';
        end

        ranking_xs = get_candidates(Model, dn, zn, dt, const, prm, config, 0);
        crit_tab = zeros(1,config.keep_xs);

        % covariance columns of the old points and of the candidates,
        % shared by the whole scoring loop
        cand = dt(ranking_xs, :);
        precomp.n_old = size(old_pts, 1);
        for m = 1:prm.M
            if precomp.n_old > 0
                kreq_old = stk_make_kreq(precomp.post_model{m}, old_pts);
                precomp.K_old{m} = stk_make_matcov(Model(m), dt, old_pts) ...
                    - precomp.lambda_mu_dt{m} * get(kreq_old, 'RS');
            else
                precomp.K_old{m} = zeros(size(dt, 1), 0);
            end
            kreq_cand = stk_make_kreq(precomp.post_model{m}, cand);
            precomp.K_cand{m} = stk_make_matcov(Model(m), dt, cand) ...
                - precomp.lambda_mu_dt{m} * get(kreq_cand, 'RS');
            [~, ~, ~, cov_joint] = stk_predict(precomp.post_model{m}, [old_pts; cand]);
            precomp.cov_joint{m} = double(cov_joint);
        end

        for r = 1:config.keep_xs
            pt = dt(ranking_xs(r),:);
            if size(unique(double([dn; old_pts; pt]), "rows"),1) == size(dn,1)+size(old_pts,1)+1
                precomp.cand_idx = r;
                crit_tab(r) = modif_mu(pt, old_pts, Model, dn, zn, dt, estim_current_temp, mu_n, IS_temp, const, prm, config, precomp);
            else
                crit_tab(r) = -inf;
            end
        end
        precomp.cand_idx = [];  % the simplex search evaluates free points

        [~, ind_max] = max(crit_tab);
        start_pt = dt(ranking_xs(ind_max),:);

        if optimize == 0
            newpt = start_pt;
        else
            options = optimset('MaxFunEvals', 50, 'MaxIter', 50, 'Display', 'off');
            crit_func = @(pt) -modif_mu(pt, old_pts, Model, dn, zn, dt, estim_current_temp, mu_n, IS_temp, const, prm, config, precomp);
            newpt = fminsearch(crit_func, double(start_pt), options);

            if ~ismember(newpt, prm.BOX) || size(unique(double([dn; old_pts; newpt]), "rows"),1) ~= size(dn,1)+size(old_pts,1)+1
                newpt = start_pt;
            end
        end

        old_pts = [old_pts; newpt];

    end

    pts_to_evaluate = old_pts;

end % function


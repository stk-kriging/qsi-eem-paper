% MODIF_MU expected modification of the estimated quantile set (EEM)
%
% Expected fraction of the population whose membership changes once the
% candidate batch [OLD_PTS; PT] is evaluated (Gauss-Hermite or QMC).

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

% precomp: optional precomputations (see get_point_EEM)

function mean_modif = modif_mu(pt, old_pts, Model, dn, zn, dt, estim_current, mu_n, IS, const, prm, config, precomp)

    % cache of the (deterministic) integration point sets
    persistent quad_cache

    dt_cand = [old_pts; pt];
    nb_cand = size(dt_cand, 1);
    nb_dt = size(dt, 1);

    duplicated_pts = size(unique(double([dn; dt_cand]), "rows"), 1) ...
        ~= size(dn,1) + nb_cand;
    if duplicated_pts
        mean_modif = 0;
        return
    end

    % posterior update mu_{n+r} = mu_n + a_n Z, Z ~ N(0, I_q), with
    % a_n = k_n(., U) V D^{-1/2} truncated at the effective rank of K_n(U)
    rank_tol = 1e-9;

    A = cell(1, prm.M);
    q = zeros(1, prm.M);

    for m = 1:prm.M
        if nargin < 13
            K_nm = cov_cond(Model(m), dn, zn(:,m), dt, dt_cand);
            [~, ~, ~, cov_cand] = stk_predict(Model(m), dn, zn(:,m), dt_cand);
        elseif isfield(precomp, 'cand_idx') && ~isempty(precomp.cand_idx)
            % columns precomputed by get_point_EEM
            r = precomp.cand_idx;
            K_nm = [precomp.K_old{m}, precomp.K_cand{m}(:, r)];
            sel = [1:precomp.n_old, precomp.n_old + r];
            cov_cand = precomp.cov_joint{m}(sel, sel);
        elseif isfield(precomp, 'K_old')
            % free point: only the last column is new
            kreq_pt = stk_make_kreq(precomp.post_model{m}, pt);
            K_nm = [precomp.K_old{m}, ...
                stk_make_matcov(Model(m), dt, pt) ...
                - precomp.lambda_mu_dt{m} * get(kreq_pt, 'RS')];
            [~, ~, ~, cov_cand] = stk_predict(precomp.post_model{m}, dt_cand);
        else
            kreq_cand = stk_make_kreq(precomp.post_model{m}, dt_cand);
            K_nm = stk_make_matcov(Model(m), dt, dt_cand) ...
                - precomp.lambda_mu_dt{m} * get(kreq_cand, 'RS');
            [~, ~, ~, cov_cand] = stk_predict(precomp.post_model{m}, dt_cand);
        end

        C = double(cov_cand);
        [V, d] = eig((C + C')/2, 'vector');
        [d, ord] = sort(d, 'descend');
        keep = (d > 0) & (d > rank_tol*d(1));
        q(m) = sum(keep);
        A{m} = K_nm * (V(:,ord(keep)) ./ sqrt(d(keep))');
    end

    % degenerate update
    Q = sum(q);
    if Q == 0
        mean_modif = 0;
        return
    end

    % weighted set of N vectors z_k in R^Q
    % legacy configs
    if ~isfield(config, 'qmc_size')
        config.qmc_size = 1024;
    end
    if ~isfield(config, 'nVar_batch')
        config.nVar_batch = 6;
    end

    % fewer 1-D nodes for batches of candidates
    if nb_cand > 1
        nv = config.nVar_batch;
    else
        nv = config.nVar;
    end

    if isempty(quad_cache) || quad_cache.qmc_size ~= config.qmc_size
        quad_cache = struct('qmc_size', config.qmc_size, 'nodes_1d', {{}}, ...
            'weights_1d', {{}}, 'tensor', {{}}, 'sobol', {{}});
    end

    if nv^Q <= 1000

        % Q-fold tensor product of the 1-D nodes
        if (numel(quad_cache.nodes_1d) < nv) || isempty(quad_cache.nodes_1d{nv})
            quantOpt.nbLevels = nv;
            quantOpt.useGaussHermite = 1;
            [nodes_1d, weights_1d] = quantization(0, 1, quantOpt);
            quad_cache.nodes_1d{nv} = nodes_1d;
            quad_cache.weights_1d{nv} = weights_1d';
        end

        if (size(quad_cache.tensor, 1) < nv) || (size(quad_cache.tensor, 2) < Q) ...
                || isempty(quad_cache.tensor{nv, Q})
            N = nv^Q;
            grids = cell(1, Q);
            [grids{1:Q}] = ndgrid(1:nv);
            I = zeros(N, Q);
            for j = 1:Q
                I(:,j) = grids{j}(:);
            end
            quad_cache.tensor{nv, Q} = struct( ...
                'w', prod(reshape(quad_cache.weights_1d{nv}(I), N, Q), 2), ...
                'Z_all', reshape(quad_cache.nodes_1d{nv}(I), N, Q)');
        end
        w = quad_cache.tensor{nv, Q}.w;
        Z_all = quad_cache.tensor{nv, Q}.Z_all;

    else

        % Sobol points mapped to N(0, I)
        if (numel(quad_cache.sobol) < Q) || isempty(quad_cache.sobol{Q})
            U01 = double(stk_sampling_sobol(config.qmc_size, Q, [], true));
            U01 = min(max(U01, 2^-53), 1 - 2^-53);  % keep norminv finite
            quad_cache.sobol{Q} = norminv(U01)';
        end
        Z_all = quad_cache.sobol{Q};
        w = ones(config.qmc_size, 1)/config.qmc_size;

    end

    N = size(Z_all, 2);

    Z = cell(1, prm.M);
    offset = 0;
    for m = 1:prm.M
        Z{m} = Z_all(offset + (1:q(m)), :);
        offset = offset + q(m);
    end

    % only the rows whose membership test can flip are evaluated
    % (Cauchy-Schwarz bound)

    n_x = nb_dt / config.pts_s;
    zmax = max(sqrt(sum(Z_all.^2, 1)));

    t0 = false(nb_dt, prm.M);
    flip_m = false(nb_dt, prm.M);
    for m = 1:prm.M  % assume outputs are independent
        mu0 = mu_n(:,m);
        if const(1,m) > const(2,m)
            t0(:,m) = (const(1,m) < mu0) | (const(2,m) > mu0);
        else
            t0(:,m) = (const(1,m) < mu0) & (const(2,m) > mu0);
        end
        delta = sqrt(sum(A{m}.^2, 2)) * zmax * (1 + 1e-12);  % inflated bound
        flip_m(:,m) = (abs(mu0 - const(1,m)) <= delta) ...
                    | (abs(mu0 - const(2,m)) <= delta);
    end
    in_C_0 = any(t0, 2);
    can_flip = any(flip_m, 2);

    cnt0 = sum(reshape(in_C_0, n_x, config.pts_s), 2);

    F = find(can_flip);
    if isempty(F)
        estim_new = repmat(cnt0/config.pts_s <= prm.alpha, 1, N);
    else
        in_C_F = false(numel(F), N);
        for m = 1:prm.M  % each output only visits its own flippable rows
            sel = flip_m(F,m);
            mu = mu_n(F(sel),m) + A{m}(F(sel),:)*Z{m};
            tm = repmat(t0(F,m), 1, N);
            if const(1,m) > const(2,m)
                tm(sel,:) = (const(1,m) < mu) | (const(2,m) > mu);
            else
                tm(sel,:) = (const(1,m) < mu) & (const(2,m) > mu);
            end
            in_C_F = in_C_F | tm;
        end
        pid = mod(F-1, n_x) + 1;  % x varies fastest in dt
        gather = sparse(pid, 1:numel(F), 1, n_x, numel(F));
        corr = gather * (double(in_C_F) - double(in_C_0(F)));
        estim_new = (cnt0 + corr)/config.pts_s <= prm.alpha;
    end

    mean_modif = 0;
    for k = 1:N
        mean_modif = mean_modif + w(k)*sum(IS .* abs(estim_current - estim_new(:,k)));
    end

end

% NEW_CONST computes the next constraint set C_k

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

function [const, thresh, mu_n] = new_const(Model, dt, dn, zn, thresh, beta, prm, config)

    is_valid = 0;

    [mu_n, var_n] = predict_mean_var(Model, dn, zn, dt, prm);

    if ~ isfield(prm, 'tol')  % interval constraint: dichotomy on the upper bound

        thresh_high = final_threshold(prm);
        thresh_low = thresh;
        thresh = (thresh_high+thresh_low)/2;
        const = threshold_to_const(prm, thresh);
        iterations = 0;
        max_iterations = 60;
        while is_valid ~= 1
            iterations = iterations + 1;

            estim = get_estim_safe(mu_n, var_n, beta, const, prm.alpha, config.pts_x, config.pts_s);
            if mean(estim) > 1.1*config.rho
                thresh_low = thresh;
                thresh = (thresh_high+thresh)/2;
            elseif mean(estim) < 0.9*config.rho
                thresh_high = thresh;
                thresh = (thresh_low+thresh)/2;
            else
                is_valid = 1;
            end

            if final_threshold(prm)/thresh < 1+10^(-12)
                thresh = final_threshold(prm);
                is_valid = 1;
            end
            if iterations >= max_iterations
                is_valid = 1;
            end
            const = threshold_to_const(prm, thresh);
        end


        %%%%%%%%%%%%%%%%%%%

    else  % band constraint: dichotomy on the tolerance
        thresh_high = thresh;
        thresh_low = final_threshold(prm);
        thresh = (thresh_high+thresh_low)/2;
        const = threshold_to_const(prm, thresh);
        iterations = 0;
        max_iterations = 60;

        while is_valid ~= 1
            iterations = iterations + 1;

            estim = get_estim_safe(mu_n, var_n, beta, const, prm.alpha, config.pts_x, config.pts_s);

            if mean(estim) > 1.1*config.rho
                thresh_high = thresh;
                thresh = (thresh_low+thresh)/2;
            elseif mean(estim) < 0.9*config.rho
                thresh_low = thresh;
                thresh = (thresh_high+thresh)/2;
            else
                is_valid = 1;
            end

            if thresh < (1 + 10^(-6))*final_threshold(prm)
                thresh = final_threshold(prm);
                is_valid = 1;
            end
            if iterations >= max_iterations
                is_valid = 1;
            end
            const = threshold_to_const(prm, thresh);
        end

    end
end % function


% INIT_SET computes the initial constraint set C_0

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

function [const, thresh, mu_n] = init_set(Model, dt, dn, zn, prm, config)

    is_valid = 0;

    mu_n = predict_mean_var(Model, dn, zn, dt, prm);

    max_iterations = 60;
    iterations = 0;

    while is_valid ~= 1
        iterations = iterations + 1;

        if ~ isfield(prm, 'tol')  % interval constraint: dichotomy on the upper bound
            thresh_high = final_threshold(prm);
            thresh_low = min(zn);
            thresh = (thresh_high+thresh_low)/2;
            const = threshold_to_const(prm, thresh);

            while is_valid ~= 1
                iterations = iterations + 1;
                estim = get_estim(mu_n, const, prm, config);

                if mean(estim) > 1.1*config.rho
                    thresh_low = thresh;
                    thresh = (thresh_high+thresh)/2;
                elseif mean(estim) < 0.9*config.rho
                    thresh_high = thresh;
                    thresh = (thresh_low+thresh)/2;
                else
                    is_valid = 1;
                end
                if iterations >= max_iterations
                    is_valid = 1;
                end
                const = threshold_to_const(prm, thresh);
            end


        %%%%%%%%%%%%%%%%%%%%%%

        else  % band constraint: dichotomy on the tolerance
            thresh_low = final_threshold(prm);
            thresh_high = 100*thresh_low;

            % enlarge the bracket if needed
            while (mean(get_estim(mu_n, threshold_to_const(prm, thresh_high), prm, config)) < 0.9*config.rho) ...
                    && (thresh_high < 10^9*thresh_low)
                thresh_high = 10*thresh_high;
            end

            thresh = (thresh_high+thresh_low)/2;
            const = threshold_to_const(prm, thresh);

            while is_valid ~= 1
                iterations = iterations + 1;

                estim = get_estim(mu_n, const, prm, config);

                if mean(estim) > 1.1*config.rho
                    thresh_high = thresh;
                    thresh = (thresh_low+thresh)/2;
                elseif mean(estim) < 0.9*config.rho
                    thresh_low = thresh;
                    thresh = (thresh_high+thresh)/2;
                else
                    is_valid = 1;
                end

                if thresh < (1 + 10^(-9))*final_threshold(prm)
                    thresh = final_threshold(prm);
                    is_valid = 1;
                end
                if iterations >= max_iterations
                    is_valid = 1;
                end
                const = threshold_to_const(prm, thresh);
            end
        end
    end
end % function


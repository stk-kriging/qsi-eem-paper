% GET_ESTIM_SAFE relaxed membership, from shifted posterior means

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

function [estim, temp_box] = get_estim_safe ...
    (mu_n, var_n, beta, const, a, xsize, usize)

    eval = zeros(xsize*usize,1);
    temp_box = const;

    if all(const(1,:) <= const(2,:))

        for m = 1:size(mu_n, 2)  % assume outputs are independent
            mu = mu_n(:,m);
            var = var_n(:,m);
            if nargout > 1
                temp_box(:,m) = const(:,m);
            end

            eval_temp = (const(2,m) >= mu+sqrt(var)*norminv(beta)) ...
                & (const(1,m) <= mu-sqrt(var)*norminv(beta));
            eval = eval | eval_temp;
        end

    elseif all(const(1,:) > const(2,:))

        for m = 1:size(mu_n, 2)  % assume outputs are independent
            mu = mu_n(:,m);
            var = var_n(:,m);
            if nargout > 1
                temp_box(:,m) = const(:,m);
            end

            eval_temp_2 = (const(2,m) > mu);
            eval_temp_2 = eval_temp_2 & (const(2,m) > mu+sqrt(var)*norminv(beta));

            eval_temp_1 = (const(1,m) < mu);
            eval_temp_1 = eval_temp_1 & (const(1,m) < mu-sqrt(var)*norminv(beta));

            eval_temp = (eval_temp_1 | eval_temp_2);
            eval = eval | eval_temp;
        end

    else

        error (['All the outputs must have the same orientation: either ' ...
            'const(1,m) <= const(2,m) for all m, or const(1,m) > const(2,m) ' ...
            'for all m.']);

    end

    eval = reshape(eval, xsize, usize);

    proba = sum(eval,2) / usize;
    estim = (proba <= a);

end % function


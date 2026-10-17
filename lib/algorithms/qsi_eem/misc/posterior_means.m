% POSTERIOR_MEANS posterior means through precomputed kriging weights
%
% Z = POSTERIOR_MEANS (MODEL, DN, ZN, D) returns the posterior means of
% the outputs at the rows of D, one column per output, for the constant
% trend models MODEL(m) fitted on (DN, ZN(:,m)).  The kriging weights
% are computed once per output and applied by blocks of rows, which is
% much cheaper than stk_predict when the variance is not needed.
%
% Z = POSTERIOR_MEANS (..., BLOCK_SIZE) sets the block size (default 1e5).

% Copyright Notice
%
%    Copyright (C) 2026 CentraleSupelec
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
%    qsi-eem-paper is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public
%    License for more details.
%
%    You should  have received a copy  of the GNU  General Public License
%    along with qsi-eem-paper.  If not, see <http://www.gnu.org/licenses/>.

function z = posterior_means (Model, dn, zn, d, block_size)

if nargin < 5
    block_size = 1e5;
end

n = size (d, 1);
z = zeros (n, size (zn, 2));

xi = double (dn);
for m = 1:size (zn, 2)
    [K, P] = stk_make_matcov (Model(m), xi);
    if size (P, 2) ~= 1 || any (P ~= P(1))
        error ('posterior_means expects a constant trend.');
    end
    zi = double (zn(:,m));
    L = chol ((K + K') / 2, 'lower');
    beta = (P' * (L' \ (L \ P))) \ (P' * (L' \ (L \ zi)));
    w{m} = L' \ (L \ (zi - P * beta));        %#ok<AGROW>
    trend(m) = P(1) * beta;                   %#ok<AGROW>
end

for i0 = 1:block_size:n
    i1 = min (i0 + block_size - 1, n);
    for m = 1:size (zn, 2)
        Kt = stk_make_matcov (Model(m), d(i0:i1, :), xi);
        z(i0:i1, m) = trend(m) + Kt * w{m};
    end
end

end % function

% PISTON_SMC cycle time of the piston simulation function
%
% U columns: M, S, V0, T0 (design) then P0, Ta (uncertain).

% Copyright Notice
%
%    Copyright (C) 2025, 2026 CentraleSupelec
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

function cycle_time = piston_smc (u)

M  = u(:, 1);  % x1
S  = u(:, 2);  % x2
V0 = u(:, 3);  % x3
T0 = u(:, 4);  % x4

P0 = u(:, 5);  % s1
Ta = u(:, 6);  % s2

K = 1000;

tmp = P0 .* V0 ./ T0 .* Ta;

A = P0.*S + 19.62*M - K*V0./S;

V = S/(2*K) .* (sqrt (A.^2 + 4*K * tmp) - A);

% cycle time in seconds
cycle_time = 2 .* pi .* (sqrt (M ./ (K + (S./V).^2 .* tmp)));

end % function

% OTL_SMC mid-point voltage of the OTL circuit, returned as 100 - Vm
%
% U columns: Rb1, Rb2, Rf, Rc1, Rc2 (design) then the gain beta (uncertain).

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

function Vm = otl_smc (u)

Rb1  = u(:, 1);  % x1
Rb2  = u(:, 2);  % x2
Rf   = u(:, 3);  % x3
Rc1  = u(:, 4);  % x4
Rc2  = u(:, 5);  % x5
beta = u(:, 6);  % s

Vb1 = 12.*Rb2 ./ (Rb1+Rb2);
term1a = (Vb1+0.74) .* beta .* (Rc2+9);
term1b = beta.*(Rc2+9) + Rf;
term1 = term1a ./ term1b;

term2a = 11.35 .* Rf;
term2b = beta.*(Rc2+9) + Rf;
term2 = term2a ./ term2b;

term3a = 0.74 .* Rf .* beta .* (Rc2+9);
term3b = (beta.*(Rc2+9)+Rf) .* Rc1;
term3 = term3a ./ term3b;

Vm = term1 + term2 + term3;
Vm = 100 - Vm;

end % function

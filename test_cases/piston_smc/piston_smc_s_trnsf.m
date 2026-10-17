% PISTON_SMC_S_TRNSF maps the uniform S sample to Beta(1/2, 1/2) margins

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

function s_trnsf = piston_smc_s_trnsf(s)

alpha = 1/2;
beta = 1/2;

s1 = s(:,1);
s2 = s(:,2);

s1 = (s1-90000)/20000;
s2 = (s2-290)/6;

s1 = 20000*betainv(double(s1), alpha, beta)+90000;
s2 = 6*betainv(double(s2), alpha, beta)+290;
s_trnsf = [s1, s2];

end % function


% PAPER_COLOR_RAMP  Light-to-dark colour ramp for ordered parameter values
%
% c = paper_color_ramp (n) returns one colour per row, derived from the
% QSI-EEM blue so that parameter sweeps share their look across figures.

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
function c = paper_color_ramp (n)

base = paper_algo_colors ('QSI_EEM');

if n == 1
    c = base;
    return
end

% 0.35 rather than a paler tint: thin lines on white need the contrast
light = base + (1 - base) * 0.35;
dark = base * 0.5;

w = linspace (0, 1, n)';
c = (1 - w) * light + w * dark;

end % function

% PAPER_ALGO_COLORS  Colour of an algorithm, shared by all the figures
%
% col = paper_algo_colors (name)

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
function col = paper_algo_colors (name)

switch char (name)
    case {'QSI_EEM', 'QSI_EEM_smc'},  col = [0.00 0.45 0.74];
    case {'QSI_m', 'QSI_SUR'},        col = [0.49 0.18 0.56];
    case {'joint_m', 'joint_SUR'},    col = [0.30 0.75 0.93];
    case 'Ranjan',                    col = [0.93 0.69 0.13];
    case 'misclassification',         col = [0.64 0.08 0.18];
    case 'random',                    col = [0.00 0.00 0.00];
    case 'BSS',                       col = [0.85 0.33 0.10];
    otherwise
        error ('Unknown algorithm: %s', char (name));
end

end % function

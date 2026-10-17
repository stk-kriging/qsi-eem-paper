% PAPER_NEW_FIGURE  Off-screen figure at the inclusion size of the paper.

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
function fig = paper_new_figure (opts)

fig = figure ('Visible', 'off', 'Color', 'white', 'Renderer', 'painters', ...
    'Units', 'centimeters', ...
    'Position', [2 2 opts.fig_width_cm opts.fig_height_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', [opts.fig_width_cm opts.fig_height_cm], ...
    'PaperPosition', [0 0 opts.fig_width_cm opts.fig_height_cm]);

% since R2025a figures follow the desktop theme, which prints dark axes
if isprop (fig, 'Theme')
    set (fig, 'Theme', 'light');
end

end % function

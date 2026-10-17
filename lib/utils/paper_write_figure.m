% PAPER_WRITE_FIGURE writes a figure as PNG (180 dpi) and EPS

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
function paper_write_figure (fig, fig_dir, fig_name)

fig_name = char (fig_name);

png_file = fullfile (fig_dir, [fig_name '.png']);
eps_file = fullfile (fig_dir, [fig_name '.eps']);

print (fig, '-dpng', '-r180', png_file);
fprintf ('Wrote %s\n', png_file);

print (fig, '-depsc2', eps_file);
fprintf ('Wrote %s\n', eps_file);

end % function


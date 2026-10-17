% PAPER_SAVE_FIGURE  Write a figure as .fig and .eps, at the size of the paper.
%
% Unlike PAPER_WRITE_FIGURE, which writes PNG + EPS, this keeps the
% reopenable .fig and a BoundingBox covering the whole canvas.

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
function paper_save_figure (fig, fig_dir, fig_name)

fig_file = fullfile (fig_dir, [fig_name '.fig']);
eps_file = fullfile (fig_dir, [fig_name '.eps']);

% openfig restores a figure saved as invisible in an unusable state
set (fig, 'Visible', 'on');

try
    savefig (fig, fig_file);
catch err
    if exist ('OCTAVE_VERSION', 'builtin') == 5
        undef = strcmp (err.identifier, 'Octave:undefined-function');
    else
        undef = strcmp (err.identifier, 'MATLAB:UndefinedFunction');
    end
    if undef
        hgsave (fig, fig_file);
    else
        rethrow (err);
    end
end

fprintf ('Wrote %s\n', fig_file);

print (fig, '-depsc2', '-painters', '-loose', eps_file);  % -loose: full canvas
fprintf ('Wrote %s\n', eps_file);

end % function

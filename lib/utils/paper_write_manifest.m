% PAPER_WRITE_MANIFEST writes manifest.txt describing the outputs of a script

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
function paper_write_manifest (out_dir, lines)

if ~ exist (out_dir, 'dir')
  mkdir (out_dir);
end

fid = fopen (fullfile (out_dir, 'manifest.txt'), 'w');
if fid < 0
  error ('Could not write manifest in %s', out_dir);
end
cleanup = onCleanup (@() fclose(fid));

fprintf (fid, 'Generated: %s\n', datestr (now (), 31));
for i = 1:(numel (lines))
  fprintf (fid, '%s\n', lines{i});
end

end % function


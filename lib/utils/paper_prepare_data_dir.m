% PAPER_PREPARE_DATA_DIR creates the output tree of an experiment

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
function paper_prepare_data_dir (data_dir)

dirs = cell(12, 1);
dirs{1} = 'doe_init';
dirs{2} = 'grid';
dirs{3} = fullfile('results', 'design');
dirs{4} = fullfile('results', 'param');
dirs{5} = fullfile('results', 'time');
dirs{6} = fullfile('results', 'const');
dirs{7} = fullfile('results', 'prop');
dirs{8} = fullfile('results', 'dev');
dirs{9} = fullfile('results', 'deviations');
dirs{10} = fullfile('results', 'graphs');
dirs{11} = fullfile('results', 'figures');
dirs{12} = fullfile('results', 'smc_logs');

for i = 1:numel(dirs)
  p = fullfile(data_dir, dirs{i});
  if ~ exist(p, 'dir')
    mkdir(p);
  end
end

end % function


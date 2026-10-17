% PAPER_GENERATE_INITIAL_DESIGNS writes the initial designs of a campaign

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
function paper_generate_initial_designs (problem, config, run_ids, data_dir)

paper_prepare_data_dir(data_dir);

[prm, ~, ~] = problem ();
dim_tot = prm.dim_x + prm.dim_s;

for i = 1:numel(run_ids)
  run_id = run_ids(i);
  file_name = sprintf('doe_init_%s_%d_init.csv', prm.name, run_id);
  file_path = fullfile(data_dir, 'doe_init', file_name);

  paper_set_seed(run_id);
  design = stk_sampling_maximinlhs(config.pts_init, dim_tot, prm.BOX);
  writematrix_(double(design), file_path);
end

end % function


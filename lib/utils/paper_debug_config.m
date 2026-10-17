% PAPER_DEBUG_CONFIG returns the configuration overrides that make an
% experiment run in seconds.

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
function config = paper_debug_config ()

config.pts_x = 40;
config.pts_s = 32;
config.nVar = 3;
config.nVar_batch = 3;
config.keep_x = 12;
config.keep_xs = 20;
config.keep = 50;
config.keep2 = 20;
config.nTraj = 10;
config.size_stop_test = 20;
config.qmc_size = 64;
config.count_resample = 2;
config.MH_steps = 2;
config.pts_x_deter = 128;
config.pts_s_deter = 32;
config.pts_eval_x = 128;
config.pts_eval_s = 32;
config.stop = 0.8;
config.stop_final = 0.8;
config.max_steps = 2;
config.T = 2;
config.axT = 1;

end % function

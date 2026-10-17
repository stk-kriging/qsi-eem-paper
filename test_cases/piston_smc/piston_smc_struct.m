% PISTON_SMC_STRUCT parameters, function and S transform of the piston case

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

function [prm, f, s_trnsf] = piston_smc_struct()

f = @piston_smc;

s_trnsf = @piston_smc_s_trnsf;

prm = struct();

prm.M = 1;

prm.dim_x = 4;
prm.dim_s = 2;

prm.BOX = stk_hrect([[30;60],[0.005;0.02], [0.002; 0.01], [340; 360], [90000; 110000], [290; 296]]);
prm.BOXx = stk_hrect([[30;60],[0.005;0.02], [0.002; 0.01], [340; 360]]);
prm.BOXs = stk_hrect([[90000; 110000], [290; 296]]);

prm.const = [-inf ; 1.148];
prm.alpha = 0.05;

prm.name = "piston_smc";
prm.name_sub = ["piston_smc"];

prm.list_cov = {             ...
    'stk_expcov_aniso';      ...
    'stk_materncov32_aniso'; ...
    'stk_materncov52_aniso'; ...
    'stk_gausscov_aniso'     };

end % function


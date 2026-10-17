% OTL_SMC_STRUCT parameters, function and S transform of the OTL case

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

function [prm, f, s_trnsf] = otl_smc_struct()

f = @otl_smc;

s_trnsf = @otl_smc_s_trnsf;

prm = struct();

prm.M = 1;

prm.dim_x = 5;
prm.dim_s = 1;

prm.BOX = stk_hrect([[50;150],[25;70], [0.5; 3], [1.2; 2.5], [0.25; 1.2], [50; 300]]);
prm.BOXx = stk_hrect([[50;150],[25;70], [0.5; 3], [1.2; 2.5], [0.25; 1.2]]);
prm.BOXs = stk_hrect(repmat([50; 300],1,prm.dim_s));

prm.const = [-inf ; 100-2.65];
prm.alpha = 0.05;

prm.name = "otl_smc";
prm.name_sub = ["otl_smc"];

prm.list_cov = {             ...
    'stk_expcov_aniso';      ...
    'stk_materncov32_aniso'; ...
    'stk_materncov52_aniso'; ...
    'stk_gausscov_aniso'     };

end % function


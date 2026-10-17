% PAPER_GET_TEST_CASE problem and configuration of a named test case

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

function [problem, config] = paper_get_test_case (case_name, debug_mode)

if nargin < 2
    debug_mode = false;
end

switch case_name

    case 'branin_mod'
        problem = @branin_mod_struct;
        config = branin_mod_config ();

    case 'double_camel'
        problem = @double_camel_struct;
        config = double_camel_config ();

    case 'hart4'
        problem = @hart4_struct;
        config = hart4_config ();

    case {'piston', 'piston_smc'}
        problem = @piston_smc_struct;
        config = piston_smc_config ();

    case {'trid', 'trid_smc'}
        problem = @trid_smc_struct;
        config = trid_smc_config ();

    case {'otl', 'otl_smc'}
        problem = @otl_smc_struct;
        config = otl_smc_config ();

    case {'rotor', 'rotor37'}
        problem = @rotor37_struct;
        config = rotor37_config ();

    otherwise
        error('Unknown case name: %s', case_name);

end % switch

if debug_mode

    debug_config = paper_debug_config ();

    names = fieldnames (debug_config);
    for i = 1:(numel (names))
        config.(names{i}) = debug_config.(names{i});
    end

end

end % function

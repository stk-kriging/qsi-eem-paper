% FINAL_THRESHOLD returns the threshold value reached at the final level
%
% SYNTAX: VALUE = final_threshold (PRM)
%
%   returns the value that the threshold of the sequence of constraint sets
%   takes when the final set PRM.const is reached: PRM.tol for a test case
%   that declares it, and the upper bound of the first output otherwise.
%
% See also: threshold_to_const, init_set, new_const, deterministic_set_const

% Copyright Notice
%
%    Copyright (C) 2026 CentraleSupelec
%
%    Author(s): Romain Ait Abdelmalek-Lomenech <romain.ait@centralesupelec.fr>

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

function value = final_threshold (prm)

if isfield (prm, 'tol')
    value = prm.tol;
else
    value = prm.const(2,1);
end

end % function

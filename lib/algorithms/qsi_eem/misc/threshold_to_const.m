% THRESHOLD_TO_CONST builds the constraint set associated with a threshold
%
% SYNTAX: CONST = threshold_to_const (PRM, VALUE)
%
%   returns the constraint set of the level whose threshold is VALUE.  Two
%   families of test cases are supported:
%
%      tolerance: the test case declares PRM.tol, and the critical region is
%         the complement of a band centred on zero.  The threshold is the
%         half-width of that band, and decreases down to PRM.tol;
%
%      interval: the critical region is the box whose bounds are the rows of
%         PRM.const.  The threshold is the upper bound of the first output,
%         and increases up to PRM.const(2,1).
%
%   VALUE = final_threshold (PRM) returns PRM.const in both cases.
%
% See also: final_threshold, init_set, new_const

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

function const = threshold_to_const (prm, value)

if isfield (prm, 'tol')

    const = (value / prm.tol) * prm.const;

else

    if any (prm.const(1,:) >= prm.const(2,:))
        error (['A critical region defined by the complement of a band ' ...
            'requires the test case to declare prm.tol.']);
    end

    const = prm.const;
    const(2,1) = value;

end

end % function

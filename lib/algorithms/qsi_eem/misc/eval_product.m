% EVAL_PRODUCT evaluates a function on a product design
%
% SYNTAX: Z = eval_product (F, XT, ST)
%
%   returns the values of F on the product design {rows of XT} x {rows of
%   ST}, ordered as ADAPT_SET orders it (the index of XT varies fastest).
%
%   A function of one argument is evaluated as F (adapt_set (XT, ST)).  A
%   test case whose evaluation can exploit the product structure declares
%   a two-argument form instead, F (XT, ST), which must return exactly the
%   same matrix; the dispatch is on nargin (F).  The shortcut being
%   carried by the function handle itself, a function that does not have
%   one -- another test case, or the posterior means of a GP model --
%   cannot end up paired with the shortcut of a different function.
%
% See also: adapt_set

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

function z = eval_product (f, xt, st)

if nargin (f) >= 2
    % the function evaluates the product design itself
    z = f (xt, st);
else
    z = f (adapt_set (xt, st));
end

end % function

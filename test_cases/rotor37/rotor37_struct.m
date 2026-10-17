% ROTOR37_STRUCT parameters, function and S transform of the ROTOR37 case
%
% See test_cases/rotor37/README.md. The outputs are relative deviations
% from the baselines, so the constraint reads directly in tolerance
% units.

function [prm, f, s_trnsf] = rotor37_struct()

f = @rotor37_qsi;

s_trnsf = @rotor37_s_trnsf;

prm = struct();

prm.M = 2;

prm.dim_x = 13;
prm.dim_s = 5;

prm.BOX = stk_hrect(prm.dim_x + prm.dim_s);
prm.BOXx = stk_hrect(prm.dim_x);
prm.BOXs = stk_hrect(prm.dim_s);

% The critical region is C = {|z1| > tol or |z2| > tol}, in relative
% deviation from the baselines.  Declaring prm.tol is what makes the SMC
% levels vary the tolerance itself, from a loose value down to prm.tol.
prm.tol = 0.18;
prm.const = prm.tol * [[1; -1], [1; -1]];

prm.alpha = 0.05;

prm.name = "rotor37";
prm.name_sub = ["rotor37"];

prm.list_cov = {             ...
    'stk_expcov_aniso';      ...
    'stk_materncov32_aniso'; ...
    'stk_materncov52_aniso'; ...
    'stk_gausscov_aniso'     };

end % function


function z = rotor37_qsi (x, s)
% MassFlow and PressureRatio, as relative deviations from their baselines.
%
% With one argument, x is a design matrix with dim_x + dim_s columns.  The
% two-argument form evaluates the product design {rows of x} x {rows of s}
% (see eval_product), exploiting the separability of the metamodel kernel:
% what depends on s only is computed once.

persistent s_cache h_cache

if nargin < 2
    z = double (rotor37 (x, 'reduced'));
    z = z(:, 1:2);
    return
end

if isempty (s_cache) || ~ isequal (s_cache, double (s))
    s_cache = double (s);
    h_cache = rotor37_prepare_s (s_cache, 'reduced');
end

zs = rotor37_eval_x (h_cache, x);

% adapt_set orders the product design with the x index varying fastest
z = [reshape(zs.MassFlow.', [], 1), reshape(zs.PressureRatio.', [], 1)];

end % function

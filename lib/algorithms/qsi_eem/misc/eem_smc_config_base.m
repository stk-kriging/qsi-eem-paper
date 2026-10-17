% EEM_SMC_CONFIG_BASE returns the usual configuration of the EEM-SMC test cases.

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
function config = eem_smc_config_base (prm)

config = struct(); %define structure

config.pts_init = 10 * (prm.dim_x + prm.dim_s); %number of points in the initial design (used by generate_doe_init.m and extract_deviation.m)

config.pts_x = 250 ; % number of points sampled in X, before importance sampling
config.pts_s = 2^9; % number of points sampled in S

config.nVar = 10; % number of points of the Gaussian quadrature used to evaluate some criteria
config.nVar_batch = 6; % same, per dimension, for batches of candidates (r > 1)

config.lognugget = log(10^-6); % log-variance of the nugget effect
config.estim_param_steps = 1; % number of steps between two estimation of the GP parameters.
                              % If == 1, parameters are estimated at every
                              % step

config.keep_x = config.pts_x; % number of points to keep in X (all of them; enforced in QSI_EEM_SMC)
config.keep_xs = 100; % number of candidate points to keep (using importance sampling)
config.count_resample = 25;
config.size_stop_test = 300; % size of the test set used by the stopping criterion

config.MH_steps = 25;
config.MH_accept = 0.35;
config.MH_delta = log(2);
% initial MH proposal: standard deviation of the uniform law on each coordinate of X
config.sigma_init = diff (double (prm.BOXx)) / sqrt (12);

config.qmc_size = 1024; % number of QMC points used by modif_mu when the tensor grid is too large

config.rho = 0.35; % target proportion of particles in each new constraint set
config.kappa = 0.1; % target relaxation of the safe estimator: safe/plain ratio close to 1+kappa

config.pts_x_deter = 1000 ; % number of points sampled in X, before importance sampling
config.pts_s_deter = 2^9;
% number of points sampled in S
config.tol_deter = 0.025;
config.stop = 0.1;
config.stop_final = 0.1;

end % function

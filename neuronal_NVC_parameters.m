function [P] = neuronal_NVC_parameters(K)
% neuronal_NVC_parameters defines structure of parameters for neuronal-NVC model, 
%                 which generates blood flow (CBF) response (see neuronal_NVC_model.m). 
%                 It applies parameter values as describe Havlicek, et al.(2015) NeuroImage  
%
% INPUT:  K - Number of cortical depths
%
% OUTPUT: P - structure with all default parameters for neuronal-NVC model
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M., Roebroeck, A., Friston, K., Gardumi, A., Ivanov, D., Uludag, K.
%           (2015) Physiologically informed dynamic cousal modeling of fMRI data, NeuroImage (122), pp. 355-372%           
% EXAMPLE:
%            K = 6;
%            P = neuronal_NVC_parameters(K)
%            disp(P);
%--------------------------------------------------------------------------
if nargin<1
  K = 6;  % default number of cortical depths
end

P.K  = K;    % Number of depths

if K<10
    P.dt = 0.01; % default integration step
elseif K<20
    P.dt = 0.005; % smaller for higher number of cortical depths
else
    P.dt = 0.001;
end

% Neuronal parameter:
%--------------------------------------------------------------------------
P.sigma   = -3;
P.mu      = 1.5;
P.lambda  = 0.2;
P.Bsigma  = [];
P.Bmu     = [];
P.Blambda = [];
P.C       = eye(K);

% NVC parameters:
% --------------------------------------------------------------------------
P.c1      = 0.6;
P.c2      = 1.5;
P.c3      = 0.6;  
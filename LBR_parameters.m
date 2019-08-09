function P = LBR_parameters(K)
% LBR_parameters defines structure of parameters for laminar BOLD response 
%                (LBR) model (see LBR_model.m). It applies parameter values
%                 as describe in Table 1 of Havlicek, M. & Uludag, K. (2019) BioRxiv  
%
% INPUT:  K - Number of cortical depths
%
% OUTPUT: P - structure with all default parameters for LBR model
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M. & Uludag, K. (2019) A dynamical model of the
%            laminar BOLD response, BioRxiv, doi: https://doi.org/10.1101/609099 
%
% EXAMPLE:
%            K = 6;
%            P = LBR_parameters(K);
%            disp(P);
%--------------------------------------------------------------------------
if nargin<1
  K = 6;  % default number of cortical depths
end

%--------------------------------------------------------------------------
P.T  = 40;   % Default time-course duration (in seconds)

P.K  = K;    % Number of depths

if K<10
    P.dt = 0.01; % default integration step
elseif K<20
    P.dt = 0.005; % smaller for higher number of cortical depths
else
    P.dt = 0.001;
end

depths = linspace(0,100,2*P.K+1); % Normalized distance to the center of individual depths (in %)
P.l    = depths(2:2:end);

% LAMINAR HEMODYNAMIC MODEL:
%--------------------------------------------------------------------------
% Baseline physiological parameters:
P.V0t   = 2.5;  % Total (regional) amount of CBV0 in the gray matter (in mL) [1-6]
P.V0t_p = 1;  % Total (regional) amount of CBV0 in the pial vein (mL) [1-6]

P.w_v = 0.5;  % CBV0 fraction of microvasculature (i.e. venules here )with respect to the total amount 
P.x_v = [];   % CBV0 fraction across depths in venules 
P.x_d = [];   % CBV0 fraction across depths in ascending veins
P.s_v = 0;    % Slope of CBV increase (decrease) in venules [0-0.3]
P.s_d = 0.3;  % Slope of CBV increase in ascending vein     [0-1.5]

P.t0v = 1;    % Transit time through microvasculature(in second)
P.E0v = 0.35; % Baseline oxygen extraction fraction in venules
P.E0d = 0.35; % Baseline oxygen extraction fraction in venules
P.E0p = 0.35; % Baseline oxygen extraction fraction in venules

% Parameters describing relative relationship between physiological variable:
% CBF-CBV coupling (steady-state)
P.alpha_v = 0.3; % For venules
P.alpha_d = 0.2; % For ascending vein
P.alpha_p = 0.1; % For pial vein

% CBF-CMRO2 coupling (steady-state)
P.n = 4;         % n-ratio   (Ref. Buxton et al. (2004) NeuroImage)

% CBF-CBV dynamic uncoupling 
P.tau_v_in = 2; % For venules - inflation 
P.tau_v_de = 2; %             - deflation

P.tau_d_in = 2; % For ascending vein - inflation 
P.tau_d_de = 2; %                    - deflation

P.tau_p_in = 2; % For pial vein - inflation 
P.tau_p_de = 2; %               - deflation

% LAMINAR BOLD SIGNAL MODEL:
%--------------------------------------------------------------------------
P.TE     = 0.028;     % echo-time (in sec)

% Hematocrit fraction
P.Hct_v  = 0.35;      % For venules, Ref. Lu et al. (2002) NeuroImage
P.Hct_d  = 0.38;      % For ascending vein
P.Hct_p  = 0.42;      % For pial vein


P.B0     = 7;                 % Magnetic field strenght (in Tesla)  
P.gyro   = 2*pi*42.6*10^6;    % Gyromagnetic constant for Hydrogen
P.suscep = 0.264*10^-6;       % Susceptibility difference between fully oxygenated and deoxygenated blood

% Water proton density:
P.rho_t  = 0.89;                % For gray matter tissue 
P.rho_v  = 0.95 - P.Hct_v*0.22; % For blood (venules) Ref. Lu et al. (2002) NeuroImage
P.rho_d  = 0.95 - P.Hct_d*0.22; % For blood (ascending vein)
P.rho_p  = 0.95 - P.Hct_p*0.22; % For blood (pial vein)
P.rho_tp = 0.95;                % For gray matter tissue % CSF   

% Relaxation rates for 7 T (in sec-1)
P.R2s_t  = 34; % For gray matter tissue
P.R2s_v  = 80; % For blood (venules)
P.R2s_d  = 85; % For blood (ascending vein)
P.R2s_p  = 90; % For blood (pial vein)


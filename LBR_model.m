function [LBR,Y,LBRpial] = LBR_model(P,cbf,cmro2)
% LBR_model calculates laminar BOLD response (LBR) as in detail described by 
%           Havlicek, M. & Uludag, K. (2019) BioRxiv  
%
% INPUTS:
%       P - structure of model parameters (see LBR_parameters.m)
%
%       cbf - matrix defining laminar cerebral blood flow (CBF) response, (time,depths). (Required);
%
%       cmro2 - matrix defining laminar changes in oxygen metabolism (CMRO2), (time,depth). 
% OUTPUTS:
%       LBR - matrix containing laminar BOLD responses in perecent signal change (time,depths)
%   
%       Y - structure with all baseline and relative physiological
%       variables underlying BOLD response
%
%       LBRpial - BOLD response of the pial vein in perecent signal change (0th depth) (time,1)
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M. & Uludag, K. (2019) A dynamical model of the
%            laminar BOLD response, BioRxiv, doi: https://doi.org/10.1101/609099 
%
% EXAMPLE:
%        For steady-state:
%               K = 6;                       % Number of depths
%               P = LBR_parameters(K);       % Get parameter structure with default values
%               cbf = ones(P.T/P.dt,K)*1.6;  % Define model input (Relative blood flow across depths)
%               P.s_d = 0.4;                 % Define the slope of increase of CBV0 in the ascening vein
%               [LBR,Y] = LBR_model(P,cbf);  % Generate the LBR
%               figure(1), 
%               plot(P.l,flipud(LBR(end,:)')); % plot the LBR profile as a function of normalized cortical depth
%               xlim([0 100]); ylim([0 6]); xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%
%        For dynamic response:
%               K = 6;  
%               P.N = neuronal_NVC_parameters(K);  % consider default parameters
%               P.N.T = 30; % (in seconds)
%               P.H = LBR_parameters(K);
%               P.H.T  = P.N.T;
%               P.H.dt = P.N.dt;
%               U.u = zeros(P.N.T/P.N.dt,K);
%               dur = 2/P.N.dt; % 2 sec stimulus
%               onset     = 2/P.N.dt;
%               offset    = onset + dur;
%               U.u(onset:offset,:) = 1;
%               [neuro, cbf]  = neuronal_NVC_model(P.N,U);
%               P.H.alpha_v   = 0.35;
%               P.H.alpha_d   = 0.2;
%               P.H.tau_d_de  = 30;
%               [LBR,Y]       = LBR_model(P.H,cbf);
%               time_axis = [0:P.H.dt:P.H.T-P.H.dt];
%               figure(1), 
%               subplot(121), plot(time_axis,cbf); xlim([time_axis(1), time_axis(end)]); ylim([0.5 2]); 
%                            xlabel('1 - Cortical depth (%)'); ylabel('Relative CBF (%)'); axis square;
%               subplot(122), plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-1 4]);  %                           xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%
%--------------------------------------------------------------------------
% 
if nargin<3
    cmro2 = [];
elseif nargin<2
    error('Not enough inputs is provided')
end


%%
% Hemodynamic model parameters:
%--------------------------------------------------------------------------
K         = P.K;      % Number of depths (Reference: Superfical depth close to CSF is k = 1;)

% BASELINE PARAMETERS:

V0t       = P.V0t;    % Total amount of CBV0 within GM tisse (in mL)
V0t_p     = P.V0t_p;  % Total amount of CBV0 in pial vein (in mL)

w_v       = P.w_v;    % Fraction of CBV0 in venules with respect to the total, 
w_d       = 1-w_v;    % Fraction of CBV0 in ascending vein with respect to the total, 

s_v       = P.s_v;    % Slope of CBV0 increase towards the surface in venules 
s_d       = P.s_d;    % Slope of CBV0 increase towards the surface in ascending vein 

% Depth-specific CBV0:
if length(P.x_v) == K,             % For venules
    x_v  = P.x_v;                  % Depth-specific fractions defined by user
else
    x_v  = 10+s_v*flipud(P.l(:));  % Possibility to define linear increase (default s_v = 0)
end
x_v      = x_v./sum(x_v);          % Fraction of CBV0 across depths in venules 

if length(P.x_v) == K,             % For ascending vein
    x_d  = P.x_d;                  % Depth-specific fractions defined by user
else
    x_d  = 10+s_d*flipud(P.l(:));  % Possibility to define linear increase 
end
x_d      = x_d./sum(x_d);          % Fraction of CBV0 across depths in venules 

V0v      = V0t*w_v*x_v;            % CBV0 in venules
V0d      = V0t*w_d*x_d;            % CBV0 in ascending vein
V0p      = V0t_p;                  % CBV0 in pial vein

% Transit time through venules (or microvasculature in general)
if length(P.t0v) == K,
    t0v    = P.t0v;              % depth-specific defined by user
else
    t0v    = ones(K,1).*P.t0v;   % default
end;

% Depth-specific baseline CBF:  
F0v     = V0v./t0v;          % Note: can be also defined directly and t0v calculated from V0v and F0v
F0d     = flipud(cumsum(flipud(F0v)));
F0p     = F0d(1);

% Depth-specific transit time:  
t0v     = V0v./F0v;   
t0d     = V0d./F0d;
t0p     = V0p./F0p;

% (check) Total mean transit time:
tt0v    = mean(t0v);   
tt0d    = mean(cumsum(t0d));
tt0     = tt0v + tt0d; % It must equal V0t./sum(F0v)

% Baseline oxygen extraction fraction
if length(P.E0v) == K,
   E0v        = P.E0v;     % depth-specific defined by user
else
   E0v        = ones(K,1).*P.E0v;     % default
end
if length(P.E0d) == K,
   E0d        = P.E0d;     % depth-specific defined by user
else
   E0d        = ones(K,1).*P.E0d;    % default  
end
E0p        = P.E0p;

%% PARAMETERS DESCRIBING RELATIVE RELATIONSHIPS BETWEEN PHYSIOLOGICAL VARIABLES:
%
% n-ratio (= (cbf-1)./(cmro2-1)). Not used if cmro2 response is directly specified as an input
if length(P.n) == K,             % For venules (microvasculature)
    n      = P.n;                % Depth-specific defined by user
else
    n      = ones(K,1)*P.n;      % Default
end;

% Grubb's exponent alpha (i.e CBF-CBV steady-state relationship)
if length(P.alpha_v) == K,       % For venules
    alpha_v    = P.alpha_v;             % Depth-specific defined by user 
else
    alpha_v    = ones(K,1).*P.alpha_v;  % Default
end;
if length(P.alpha_d) == K,       % For ascending vein
    alpha_d    = P.alpha_d;             % Depth-specific defined by user  
else
    alpha_d    = ones(K,1).*P.alpha_d;  % Default
end
alpha_p        = P.alpha_p;      % For pial vein

% CBF-CBV uncoupling (tau) during inflation and deflation:
if length(P.tau_v_in) == K,      % For venules (inflation)
    tau_v_in  = P.tau_v_in;             % Depth-specific defined by user
else
    tau_v_in  = ones(K,1).*P.tau_v_in;  % Default  
end
if length(P.tau_v_de) == K,      % For  venules (deflation)
    tau_v_de  = P.tau_v_de;             % Depth-specific defined by user  
else
    tau_v_de  = ones(K,1).*P.tau_v_de;  % Default  
end
if length(P.tau_d_in) == K,      % For ascending vein (inflation)
    tau_d_in  = P.tau_d_in;             % Depth-specific defined by user 
else
    tau_d_in  = ones(K,1)*P.tau_d_in;   % Default  
end;
if length(P.tau_d_de) == K,       % For ascending vein (deflation)
    tau_d_de  = P.tau_d_de;             % Depth-specific defined by user 
else
    tau_d_de  = ones(K,1).*P.tau_d_de;  % Default
end
tau_p_in      = P.tau_p_in;       % For pial vein (inflation)
tau_p_de      = P.tau_p_de;       % For pial vein (deflation)


%% Parameters for laminar BOLD signal equation (for 7 T field strenght):
%--------------------------------------------------------------------------
% Baseline CBV in fraction with respect to GM tissue
V0vq = V0v./100*K;
V0dq = V0d./100*K;
V0pq = V0p./100*K;

TE     = P.TE;          % echo-time (sec) 

Hct_v  = P.Hct_v;       % Hematocrit fraction
Hct_d  = P.Hct_d;
Hct_p  = P.Hct_p;
B0     = P.B0;          % Field strenght        
gyro   = P.gyro;        % Gyromagnetic constant 
suscep = P.suscep;      % Susceptibility difference

nu0v   = suscep*gyro*Hct_v*B0;
nu0d   = suscep*gyro*Hct_d*B0;
nu0p   = suscep*gyro*Hct_p*B0; 

% Water proton density 
rho_t  = P.rho_t;  % In GM tissue
rho_v  = P.rho_v;  % In blood (venules) Ref. Lu et al. (2002) NeuroImage
rho_d  = P.rho_d;  % In blood (ascening vein) 
rho_p  = P.rho_p;  % In blood (pial vein) 
rho_tp = P.rho_tp; % In in tissue and CSF 

% Relaxation rates (in sec-1):
if length(P.R2s_t) == K,  % For tissue
    R2s_t  = P.R2s_t;           %
else
    R2s_t  = ones(K,1).*P.R2s_t;   % (sec-1)
end
if length(P.R2s_v) == K,  % For venules
    R2s_v  = P.R2s_v;               % (sec-1)
else
    R2s_v  = ones(K,1)*P.R2s_v; % (sec-1) 
end
if length(P.R2s_d) == K,  % For ascening vein
    R2s_d  = P.R2s_d;           % (sec-1)
else
    R2s_d  = ones(K,1)*P.R2s_d; % (sec-1)  
end
R2s_p  = P.R2s_p;         % For pial vein 

% (Baseline) Intra-to-extra-vascular signal ratio
ep_v   = rho_v./rho_t.*exp(-TE*R2s_v)./exp(-TE*R2s_t);        % For venules
ep_d   = rho_d./rho_t.*exp(-TE*R2s_d)./exp(-TE*R2s_t);        % For ascending vein
ep_p   = rho_p./rho_tp.*exp(-TE*R2s_p)./exp(-TE*R2s_t);       % For pial vein 

% Slope of change in R2* of blood with change in extraction fration during activation 
r0v    = 228;      % For venules   
r0d    = 232;      % For ascending vein
r0p    = 236;      % For pial vein

H0     = 1./(1 - V0vq - V0dq + ep_v.*V0vq + ep_d.*V0dq);  % constant in front
H0p    = 1./(1 - V0pq + ep_p.*V0pq);

k1v     = 4.3.*nu0v.*E0v.*TE;
k2v     = ep_v.*r0v.*E0v.*TE;
k3v     = 1 - ep_v;

k1d     = 4.3.*nu0d.*E0d.*TE;
k2d     = ep_v.*r0d.*E0d.*TE;
k3d     = 1 - ep_d;

k1p     = 4.3.*nu0p.*E0p.*TE;
k2p     = ep_p.*r0p.*E0p.*TE;
k3p     = 1 - ep_p;


%% Initial conditions:

Xk       = zeros(K,4);
Xp       = zeros(1,2);

yk       = Xk;
yp       = Xp;

f_d      = ones(K,1);
dv_d     = ones(K,1);
dHb_d    = ones(K,1);

tau_v    = tau_v_in;
tau_d    = tau_d_in;
tau_p    = tau_p_in;

% integration step
dt = P.dt;

LBR       = zeros(P.T/dt,K);
LBRpial   = zeros(P.T/dt,K);
%%
for t = 1:P.T/dt

    Xk      = exp(Xk);    % log-normal transformation (Stephan et al.(2008), NeuroImage)
    Xp      = exp(Xp);
    
    % model input (laminar CBF response):
    f_a = cbf(t,:)';
    
    % VENULES COMPARTMENTS:
    %--------------------------------------------------------------------------
    % blood outflow from venules compartment
    if sum(alpha_v)>0
        f_v     = (V0v.*Xk(:,1).^(1./alpha_v) + F0v.*tau_v.*f_a)./(V0v+F0v.*tau_v);
    else
        f_v     = f_a;
    end
    % change in blood volume in venules:
    dv_v        = (f_a - f_v)./t0v;
    % change in oxygen matabolims (CMRO2)
    if isempty(cmro2)
        m        = (f_a + n-1)./n;  %(if not specified directly)
    else
        m        = cmro2(t,:)';
    end
    % change in deoxyhemoglobin content venules:
    dHb_v        = (m - f_v.*Xk(:,2)./Xk(:,1))./t0v;


    % ASCENDING VEIN COMPARTMENTS:
    %--------------------------------------------------------------------------    
    % blood outflow from Kth depth of ascending vein compartment (deepest depth):
    if alpha_d(end)>0
        f_d(end)  = (V0d(end).*Xk(end,3).^(1./alpha_d(end)) + tau_d(end).*f_v(end).*F0v(end))./(V0d(end)+F0d(end).*tau_d(end));
    else
        f_d(end)  = f_v(end)*F0v(end)./F0d(end);
    end
    % changes in blood volume and deoxyhemoglobin in ascending vein (deepest depth):
    dv_d(end)     = (f_v(end) - f_d(end))./t0d(end);
    dHb_d(end)    = (f_v(end).*Xk(end,2)./Xk(end,1) - f_d(end).*Xk(end,4)./Xk(end,3))./t0d(end);
    
    % blood outflow from other comparments of ascending vein:
    for i = K-1:-1:1,
        if alpha_d(i)>0
            f_d(i)     = (V0d(i).*Xk(i,3).^(1./alpha_d(i)) + tau_d(i).*(f_v(i).*F0v(i)+f_d(i+1).*F0d(i+1)))./(V0d(i)+F0d(i).*tau_d(i));
        else
            f_d(i)     = f_v(i)*F0v(i)./F0d(i)+f_d(i+1)*F0d(i+1)./F0d(i);
        end
        % changes in blood volume and deoxyhemoglobin in ascending vein:
        dv_d(i)    = (f_v(i).*F0v(i)./F0d(i) + f_d(i+1).*F0d(i+1)./F0d(i) - f_d(i))./t0d(i);
        dHb_d(i)   = (f_v(i).*F0v(i)./F0d(i).*Xk(i,2)./Xk(i,1) + f_d(i+1).*F0d(i+1)./F0d(i).*Xk(i+1,4)./Xk(i+1,3) - f_d(i).*Xk(i,4)./Xk(i,3))./t0d(i);

    end;
    
    % PIAL VEIN COMPARTMENT:
    %--------------------------------------------------------------------------    

    % blood outflow from pial vein:
    if alpha_p>0
        f_p     = (V0p.*Xp(1).^(1./alpha_p) + F0p.*tau_p.*f_d(1))./(V0p+F0p.*tau_p);
    else
        f_p     = f_d(1);
    end;
    % changes in blood volume and deoxyhemoglobin in pial vein:
    dv_p  = (f_d(1) - f_p)./t0p;
    dHb_p = (f_d(1).*Xk(1,4)./Xk(1,3) - f_p.*Xp(2)./Xp(1))./t0p;
    
    
    % Intergrated changes to previous time point
    yk(:,1)  = yk(:,1) + dt*(dv_v./Xk(:,1));
    yk(:,2)  = yk(:,2) + dt*(dHb_v./Xk(:,2));
    yk(:,3)  = yk(:,3) + dt*(dv_d./Xk(:,3));
    yk(:,4)  = yk(:,4) + dt*(dHb_d./Xk(:,4));
    
    yp(:,1)  = yp(:,1) + dt*(dv_p./Xp(1));
    yp(:,2)  = yp(:,2) + dt*(dHb_p./Xp(2));
    


    Xk        = yk;
    Xp        = yp;
    
    tau_v     = tau_v_in;
    tau_d     = tau_d_in;
    tau_p     = tau_p_in;
 
    % check for deflation (negative derivative)
    tau_v(dv_v<0)  = tau_v_de(dv_v<0);
    tau_d(dv_d<0)  = tau_d_de(dv_d<0);
    tau_p(dv_p<0)  = tau_p_de(dv_p<0);
    
    
    % venules:
    m_v  = m;
    v_v  = exp(yk(:,1)); % log-normal transformation
    q_v  = exp(yk(:,2));
    % draining vein:
    v_d  = exp(yk(:,3));
    q_d  = exp(yk(:,4));
    % pail vein:
    v_p  = exp(yp(:,1));
    q_p  = exp(yp(:,2));
    
    
    % save physiological variable:
    Y.fa(t,:) = f_a;
    Y.mv(t,:) = m_v;
    Y.qv(t,:) = q_v;
    Y.qd(t,:) = q_d;
    Y.qp(t,:) = q_p;

    Y.vv(t,:) = v_v;
    Y.vd(t,:) = v_d;
    Y.vp(t,:) = v_p;

    
    
    LBR(t,:) = H0.*((1-V0vq-V0dq).*(k1v.*V0vq.*(1-q_v)       +k1d.*V0dq.*(1-q_d))+...
                                    +k2v.*V0vq.*(1-q_v./v_v) +k2d.*V0dq.*(1-q_d./v_d)+...
                                    +k3v.*V0vq.*(1-v_v)      +k3d.*V0dq.*(1-v_d)).*100;
    
    
    LBRpial(t,:) = H0p.*((1-V0pq).*(k1p.*V0pq.*(1-q_p))      +k2p.*V0pq.*(1-q_p./v_p)+...
                                    k3p.*V0pq.*(1-v_p)).*100;
                         
end;

% save baseline physiological parameters
Y.F0v  = F0v;
Y.F0d  = F0d;
Y.F0p  = F0p;

Y.V0v  = V0v;
Y.V0d  = V0d;
Y.V0p  = V0p;

Y.V0vq = V0vq;
Y.V0dq = V0dq;
Y.V0pq = V0pq;

Y.t0v  = t0v;
Y.t0d  = t0d;
Y.t0p  = t0p;
Y.tt0v = tt0v;
Y.tt0d = tt0d;
Y.tt0  = tt0;

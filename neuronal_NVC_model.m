function [neuro,cbf] = neuronal_NVC_model(P,U),
% neuronal_NVC_model is a simple model to generate neuronal response (as balance 
%                 between excitatory and inhibitory activity) followed by blood flow (CBF)
%                 response, based on Havlicek, et al.(2015) NeuroImage  
%
% INPUT:  K - Number of cortical depths
%
% OUTPUT: Y - structure with all baseline and relative physiological
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M., Roebroeck, A., Friston, K., Gardumi, A., Ivanov, D., Uludag, K.
%           (2015) Physiologically informed dynamic cousal modeling of fMRI data, NeuroImage (122), pp. 355-372
%
% EXAMPLE:
%            K = 6;
%            P = neuronal_NVC_parameters(K)
%            disp(P);
%--------------------------------------------------------------------------
K      = P.K;

% Neuronal parameter:
%--------------------------------------------------------------------------
sigma   = P.sigma;     % self-inhibitory connection 
mu      = P.mu;        % inhibitory-excitatory connection
lambda  = P.lambda;    % inhibitory gain
Bsigma  = P.Bsigma;    % modulatory parameter of self-inhibitory connection
Bmu     = P.Bmu;       % modulatory parameter of inhibitory-excitatory connection
Blambda = P.Blambda;   % modulatory parameter of inhibitory connection
C       = P.C;
% NVC parameters:
% --------------------------------------------------------------------------
c1      = P.c1;
c2      = P.c2;
c3      = P.c3;

% Initial condtions:
Xn  = zeros(K,4);
yn  = Xn;

dt  = P.dt;
neuro = zeros(P.T/dt,K);
cbf   = zeros(P.T/dt,K);
for t = 1:P.T/dt
    Xn(:,4) = exp(Xn(:,4));

    A = eye(K)*sigma;
    MU = ones(K,1)*mu;
    LAM = ones(K,1)*lambda;
    for i = 1:size(Bsigma,2)
        A = A + diag(Bsigma(:,i)).*U.m(t,i);
    end;
    for i = 1:size(Bmu,2)
        MU = MU + Bmu(:,i).*U.m(t,i);
    end;
    for i = 1:size(Blambda,2)
        LAM = LAM + Blambda(:,i).*U.m(t,i);
    end
    %----------------------------------------------------------------------
    % Neuronal (excitatiry & inhibitory)
    yn(:,1)   = yn(:,1) + dt*(A*Xn(:,1) - MU.*Xn(:,2) + C*U.u(t,:)');

    yn(:,2)   = yn(:,2) + dt*(LAM.*(-Xn(:,2) +  Xn(:,1)));
    %----------------------------------------------------------------------
    % Vasoactive signal:
    yn(:,3)   = yn(:,3) + dt*(Xn(:,1) - c1.*(Xn(:,3)));
    %----------------------------------------------------------------------
    % Inflow:
    df_a      = c2.*Xn(:,3) - c3.*(Xn(:,4)-1);
    yn(:,4)   = yn(:,4) + dt*(df_a./Xn(:,4));
    
    Xn         = yn;
     
    cbf(t,:)   = exp(yn(:,4))';
    neuro(t,:) = yn(:,1)';
end
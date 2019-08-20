% Example demonstrating dynamic laminar BOLD response (LBR)for short (2 sec) and longer
% (20 sec) stimulus durations.

close all; clear all;

set(0,'DefaultAxesFontSize', 14, ...
      'defaultLineLineWidth', 2, ...
      'defaultLineMarkerSize',15,...
      'DefaultAxesTitleFontWeight', 'normal');      
  
K         = 6;   % number of depths

% First example: Laminar BOLD response to short 2 sec stimulus
%==========================================================================
% Specify neuronal and NVC model:
%--------------------------------------------------------------------------
P.N       = neuronal_NVC_parameters(K);  % get default parameters (see inside the function)
P.N.T     = 30;               % Total lenght of the response (in seconds)
dur       = 2/P.N.dt;         % Stimulus duration (in second, e.g. 2 sec) ... dt - refers to integration step
onset     = 3/P.N.dt;         % Stimulus onset time (in seconds) 
offset    = onset + dur;      % Stimulus offset time (in seconds) 
U.u       = zeros(P.N.T/P.N.dt,K);   % Matrix with input vectors to the neuronal model (one column per depth)
U.u(onset:offset,:) = 1;             % Set one during stimulus window
[neuro, cbf]  = neuronal_NVC_model(P.N,U);  % Generate the neuronal and cerebral blood flow response (CBF)
  
% Specify LBR model:
%--------------------------------------------------------------------------  
P.H       = LBR_parameters(K); % get default parameters (see inside the function), 
                               % NOTE: By default baseline CBV is increasing towards the surface in the ascending vein
P.H.T     = P.N.T;             % copy the lenght of the response from neuronal specification
  
P.H.alpha_v   = 0.35;          % Choose steady-state CBF-CBV coupling for venules
P.H.alpha_d   = 0.2;           % Choose steady-state CBF-CBV coupling for ascending vein
P.H.tau_d_de  = 30;            % Choose dynamic CBF-CBV uncoupling for ascending vein

[LBR,Y] = LBR_model(P.H,cbf);  % Generate the laminar bold response


time_axis = [0:P.H.dt:P.H.T-P.H.dt] - onset*P.N.dt; % time axis in seconds

% Display underlying physiological responses
figure(1),
subplot(231), plot(time_axis,cbf); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CBF in MV (%)'); axis square; 
subplot(232), plot(time_axis,Y.mv); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CMRO_2 in MV (%)'); axis square;
subplot(233), plot(time_axis,Y.vv); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CBV in MV (%)'); axis square;
subplot(234), plot(time_axis,Y.qv); xlim([time_axis(1), time_axis(end)]); ylim([0.7 1.2]);
xlabel('Time (s)'); ylabel('Relative dHb in MV (%)'); axis square;
subplot(235), plot(time_axis,Y.vd); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CBV in AV (%)'); axis square;
subplot(236), p = plot(time_axis,Y.qd); xlim([time_axis(1), time_axis(end)]); ylim([0.7 1.2]);
xlabel('Time (s)'); ylabel('Relative dHb in AV (%)'); axis square; legend([p(1) p(end)],{'Upper','Lower'});
 

% Display laminar BOLD response
figure(2),
subplot(131), p = plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-1 4]);  %                         
xlabel('Time (s)'); ylabel('LBR (%)'); axis square;  legend([p(1) p(end)],{'Upper','Lower'});

% calculate time to peak (TTP) and time to undershoot (TTU) with respect to
% the stimulus onset and offset, respectively
[Peak_Amp,Peak_Pos] = max(LBR(onset:end,:));
[PSU_Amp,PSU_Pos]   = min(LBR(offset:end,:));
TTP = time_axis(onset+Peak_Pos)';
TTU = time_axis(offset+PSU_Pos)'-(offset-onset)*P.N.dt;  
% Display TTP and TTU as function of cortical depth
subplot(132), plot(P.H.l,flipud(TTP),'.-'); xlim([0 100]); ylim([0 12]);  %                          
xlabel('1 - Cortical depth (%)'); ylabel('TTP (s)'); axis square;
subplot(133), plot(P.H.l,flipud(TTU),'.-'); xlim([0 100]); ylim([0 12]);  %                         
xlabel('1 - Cortical depth (%)'); ylabel('TTU (%)'); axis square;
  



% Second example: Laminar BOLD response to short 20 sec stimulus
%==========================================================================
% Specify neuronal and NVC model:
%--------------------------------------------------------------------------
P.N       = neuronal_NVC_parameters(K);  % get default parameters (see inside the function)
P.N.mu    = 0.5;                % Inhibitory-Excitatory coupling (controls strenght of neuronal (CBF) response transients)
P.N.C     = 0.7;              % Stimulus strenght
P.N.T     = 60;               % Total lenght of the response (in seconds)
dur       = 20/P.N.dt;        % Stimulus duration (in second, e.g. 2 sec) ... dt - refers to integration step
onset     = 3/P.N.dt;         % Stimulus onset time (in seconds) 
offset    = onset + dur;      % Stimulus offset time (in seconds) 
U.u       = zeros(P.N.T/P.N.dt,K);   % Matrix with input vectors to the neuronal model (one column per depth)
U.u(onset:offset,:) = 1;             % Set one during stimulus window
[neuro, cbf]  = neuronal_NVC_model(P.N,U);  % Generate the neuronal and cerebral blood flow response (CBF)
  
% Specify LBR model:
%--------------------------------------------------------------------------  
P.H       = LBR_parameters(K); % get default parameters (see inside the function)
P.H.T     = P.N.T;             % copy the lenght of the response from neuronal specification
  
P.H.alpha_v   = 0.35;          % Choose steady-state CBF-CBV coupling for venules
P.H.alpha_d   = 0.2;           % Choose steady-state CBF-CBV coupling for ascending vein

P.H.tau_v_in  = 10;            % Choose dynamic CBF-CBV uncoupling for ascending vein (during inflation)
P.H.tau_d_in  = 20;            % Choose dynamic CBF-CBV uncoupling for ascending vein (during inflation)
                               % NOTE: Depending how strong are neruonal (CBF) response transients - CBF-CBV uncoupling (during inflation)in
                               % only in the ascending vein can make the TTP peak even longer in lower depths compared to the upper depths
P.H.tau_v_de  = 20;            % Choose dynamic CBF-CBV uncoupling for ascending vein (during deflation) - longer for longer stimulus durations
                  
P.H.tau_d_de  = 50;            % Choose dynamic CBF-CBV uncoupling for ascending vein (during deflation) - longer for longer stimulus durations

[LBR,Y] = LBR_model(P.H,cbf);  % Generate the laminar bold response


time_axis = [0:P.H.dt:P.H.T-P.H.dt] - onset*P.N.dt; % time axis in seconds

% Display underlying physiological responses
figure(3),
subplot(231), plot(time_axis,cbf); xlim([time_axis(1), time_axis(end)]); ylim([0.8 2]);
xlabel('Time (s)'); ylabel('Relative CBF in MV (%)'); axis square; 
subplot(232), plot(time_axis,Y.mv); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CMRO_2 in MV (%)'); axis square;
subplot(233), plot(time_axis,Y.vv); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CBV in MV (%)'); axis square;
subplot(234), plot(time_axis,Y.qv); xlim([time_axis(1), time_axis(end)]); ylim([0.7 1.2]);
xlabel('Time (s)'); ylabel('Relative dHb in MV (%)'); axis square;
subplot(235), plot(time_axis,Y.vd); xlim([time_axis(1), time_axis(end)]); ylim([0.8 1.6]);
xlabel('Time (s)'); ylabel('Relative CBV in AV (%)'); axis square;
subplot(236), p = plot(time_axis,Y.qd); xlim([time_axis(1), time_axis(end)]); ylim([0.7 1.2]);
xlabel('Time (s)'); ylabel('Relative dHb in AV (%)'); axis square; legend([p(1) p(end)],{'Upper','Lower'});
 

% Display laminar BOLD response
figure(4),
subplot(131), p = plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-2 7]);  %                         
xlabel('Time (s)'); ylabel('LBR (%)'); axis square;  legend([p(1) p(end)],{'Upper','Lower'});

% calculate time to peak (TTP) and time to undershoot (TTU) with respect to
% the stimulus onset and offset, respectively
[Peak_Amp,Peak_Pos] = max(LBR(onset:end,:));
[PSU_Amp,PSU_Pos]   = min(LBR(offset:end,:));
TTP = time_axis(onset+Peak_Pos)';
TTU = time_axis(offset+PSU_Pos)'-(offset-onset)*P.N.dt;  
% Display TTP and TTU as function of cortical depth
subplot(132), plot(P.H.l,flipud(TTP),'.-'); xlim([0 100]); ylim([0 12]);  %                          
xlabel('1 - Cortical depth (%)'); ylabel('TTP (s)'); axis square;
subplot(133), plot(P.H.l,flipud(TTU),'.-'); xlim([0 100]); ylim([0 12]);  %                         
xlabel('1 - Cortical depth (%)'); ylabel('TTU (%)'); axis square;
  
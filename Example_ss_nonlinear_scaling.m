% Steady-state example demonstrating (non)linear scaling of depth-dependent BOLD
% response with 'stimulus strenght' or incresed relative CBF;
close all; clear all;

set(0,'DefaultAxesFontSize', 14, ...
      'defaultLineLineWidth', 2, ...
      'defaultLineMarkerSize',15,...
      'DefaultAxesTitleFontWeight', 'normal');      

K = 6;                       % We will consider six depths
P{1} = LBR_parameters(K);    % Get parameter structure for LBR model
                             % By default we consider 40 sec stimulus
                             % duration in order to reach steady-state (i.e P{1}.T = 40)
% Define laminar profile relative CBF (considering six depth) with 60% at 
% for different scaling (1.2 1.4 1.6 1.6) of relarive CBF
cbf_original   = [1 1 1.1 1.1 1 1];  %!!! variable CBF across depths
% the top and low depths and 30% in the middle depths 
cbf{1}      = kron(cbf_original*1.2,ones(P{1}.T/P{1}.dt,1));
P{1}.V0t  = 2;    % define total amount of baseline CBV in the GM (in mL)
P{1}.w_v  = 0.4;  % define fraction of microvasculature (venules) with respect to the AV
P{1}.s_d  = 0.3;  % increase of baseline CBV in the AV towards the surface
 
% Call the LBR model to generate the laminar BOLD profile:
[LBR{1},Y{1}] = LBR_model(P{1},cbf{1});

% Repeat the same simulations but with different (larger) increase of baseline CBV
% in the AV towards the surface:
cbf{2}          = kron(cbf_original*1.4,ones(P{1}.T/P{1}.dt,1));
P{2}          = P{1};
[LBR{2},Y{2}] = LBR_model(P{2},cbf{2});


% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
P{3}          = P{1};
cbf{3}          = kron(cbf_original*1.6,ones(P{1}.T/P{1}.dt,1));
[LBR{3},Y{3}] = LBR_model(P{3},cbf{3});


% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
P{4}          = P{1};
cbf{4}          = kron(cbf_original*1.8,ones(P{1}.T/P{1}.dt,1));
[LBR{4},Y{4}] = LBR_model(P{4},cbf{4});

% Display results:
figure(1)

for i = 1:length(LBR)
    subplot(131),
    plot(P{1}.l,flipud(cbf{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([1 2]); 
    xlabel('1 - Cortical depth (%)'); ylabel('relative CBF (-)'); axis square; title('Laminar CBF profile')
    subplot(132),
    plot(P{i}.l,flipud(Y{i}.V0dq*100),'.-'); hold on; xlim([0 100]); ylim([0 3]); 
    xlabel('1 - Cortical depth (%)'); ylabel('Baseline CBV (%)'); axis square; title('Laminar baseline CBV profile')
    subplot(133),
    plot(P{i}.l,flipud(LBR{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([0 5]); title('Laminar BOLD profile')
    xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
    scale(i) = LBR{i}(end,1)./LBR{i}(end,end);
end;
hold off;


figure(2)
bar(scale), title('Ratio between upper and lower depth'); % if flat across different stimulus strengths
ylim([0 3])% it suggests linear. In the case, of variable relative CBF (neuronal activity) across depth it is not linear.
% But of course if the variation is very small, this nonlinearity will be negligible. 

figure(3) % relation ship between CBF and BOLD is still nonlinear
          % we also know that the relatinship between CBF and neuronal
          % response is often nonlinear (especially for more sustained responses)
          % combination of these two nonlinearities often leads to
          % observation that relationship between neuronal and BOLD
          % response is more or less linear.
plot([cbf{1}(end,1),cbf{2}(end,1),cbf{3}(end,1),cbf{4}(end,1)]',...
     [LBR{1}(end,1),LBR{2}(end,1),LBR{3}(end,1),LBR{4}(end,1)]); xlim([1 2]); hold on;
plot([cbf{1}(end,end),cbf{2}(end,end),cbf{3}(end,end),cbf{4}(end,end)]',...
      [LBR{1}(end,end),LBR{2}(end,end),LBR{3}(end,end),LBR{4}(end,end)]); xlim([1 2]); hold off
 xlabel('realative CBF (-)'), ylabel('BOLD (%)'); title('CBF vs BOLD'); legend({'Upper','Lower'})
 
 
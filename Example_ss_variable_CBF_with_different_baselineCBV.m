% Steady-state example demonstrating variable relative CBF across depths 
% with different amounts of baseline CBV increase of asceding veins (AV)
% towards the gray matter (GM) surface. 

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
% the top and low depths and 30% in the middle depths 
cbf      = kron([1.6,1.6,1.3,1.3,1.6,1.6],ones(P{1}.T/P{1}.dt,1));
P{1}.V0t = 2;    % define total amount of baseline CBV in the GM (in mL)
P{1}.w_v = 0.4;  % define fraction of microvasculature (venules) with respect to the AV
P{1}.s_d = 0.3;  % increase of baseline CBV in the AV towards the surface
 
% Call the LBR model to generate the laminar BOLD profile:
[LBR{1},Y{1}] = LBR_model(P{1},cbf);

% Repeat the same simulations but with different (larger) increase of baseline CBV
% in the AV towards the surface:
P{2}     = P{1};
P{2}.s_d = 1.1;
[LBR{2},Y{2}] = LBR_model(P{2},cbf);

% Repeat the same simulations but with a constant baseline CBV across
% depths in the AV:
P{3}     = P{1};
P{3}.s_d = 0;  % slope equals zero
[LBR{3},Y{3}] = LBR_model(P{3},cbf);

% Display results:
figure(1)
subplot(131),
    plot(P{1}.l,flipud(cbf(end,:)'),'.-'); hold on; xlim([0 100]); ylim([1 2]); 
    xlabel('1 - Cortical depth (%)'); ylabel('relative CBF (-)'); axis square; title('Laminar CBF profile')
    

for i = 1:length(LBR)
    subplot(132),
    plot(P{i}.l,flipud(Y{i}.V0dq*100),'.-'); hold on; xlim([0 100]); ylim([0 3]); 
    xlabel('1 - Cortical depth (%)'); ylabel('Baseline CBV (%)'); axis square; title('Laminar baseline CBV profile')
    subplot(133),
    plot(P{i}.l,flipud(LBR{i}(end,:)'),'.-'); hold on; xlim([0 100]); ylim([0 5]); title('Laminar BOLD profile')
    xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
end;
hold off;
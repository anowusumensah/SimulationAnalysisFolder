close all; clear; clc;

lw = 5;
stepVal = 286;
step = stepVal + 5;
xlineV = stepVal; 
figure('Color','w')
Iso_285 = textread("LQT2_Drug_286");

timeAxis  = Iso_285(:,1);
myo285 = Iso_285(:,3);
purk285 = Iso_285(:,3);


% % 

plot(Iso_285(:,1), Iso_285(:,3), 'b','LineWidth',lw); hold on
plot(Iso_285(:,1), Iso_285(:,2), 'r','LineWidth',lw);
% % xline(285, '--k', LineWidth=2)

% % 
% % Iso_287 = textread("LQT2_Iso_285");
% % plot(Iso_287(:,1), Iso_287(:,3), 'b','LineWidth',lw); hold on
% % plot(timeAxis(timeAxis < step), myo285(timeAxis < step), 'b','LineWidth',lw);
% % 
% % plot(Iso_287(:,1), Iso_287(:,2), 'r','LineWidth',lw);
% % plot(timeAxis(timeAxis < step), purk285(timeAxis < step), 'r','LineWidth',lw);
% % xline (stepVal, '--k', LineWidth=2)


xlim([min(Iso_285(:,1)), 1038])
ylim([-90  60])

figure('Color','w')
plot(Iso_285(:,1), Iso_285(:,3), 'b','LineWidth',lw); hold on
plot(Iso_285(:,1), Iso_285(:,2), 'r','LineWidth',lw);
xlim([stepVal- 5, stepVal + 20])
ylim([-90  60])


% % legend("M", "P")

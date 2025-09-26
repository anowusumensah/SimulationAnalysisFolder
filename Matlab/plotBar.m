close all; clear; clc;


% Sample data
categories = {'Ctrl', 'Ctrl_Drug', 'LQT2', 'LQT2_Drug'};

values(:, 1) = [231, 251, 251, 286];

values(:, 2) = [288.74089, 377.05987, 316.874804, 425.082804];

diffVals = values(:,2) - values(:,1); 


% Create bar chart
wz = 1;
figure('Color','w');
% % b = bar(values, 'FaceColor', 'flat');  % Use 'flat' to enable individual bar coloring
b = bar(diffVals, 'FaceColor', 'flat','BarWidth', wz);  % Use 'flat' to enable individual bar coloring
% Assign individual colors to each bar
% % b.CData(1,:) = [1 0 0];    % Red
b.CData(1,:) = [0.68, 0.85, 0.90];    % Blue
b.CData(2,:) = [0 0 1];    % Blue
b.CData(3,:) = [1.00, 0.71, 0.76];
b.CData(4,:) = [1 0 0];

% % b.CData(3,:) = [0 1 0];    % Green
% % b.CData(4,:) = [1 1 0];    % Yellow
ax = gca;
ax.TickLabelInterpreter = 'none';
ax.FontSize  = 20;
ax.FontName = 'Arial';
ax.FontWeight =  'bold';
ax.LineWidth = 2; 
% % ax.TickLength = [0.02 0.02];




% Customize chart
set(gca, 'XTickLabel', categories);  % Set x-axis labels
% % xlabel('Experiment');
% % ylabel('Conduction Delay (ms)');
box off


% % title('Fruit Sales in June');
% % grid on;

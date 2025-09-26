% File: plot_with_inset.m
% Description: Plot time vs two voltage traces with inset highlighting a time range
close all; clear; clc;

figure('Color','w')

% Load data
stemTime = 250;
data = load('LQT2_Drug_savF_1');  % Replace with your actual filename
time = data(:, 1);
voltage1 = data(:, 2);
voltage2 = data(:, 3);

% Main figure
% % figure;
plot(time, voltage1, 'r', 'LineWidth', 1.5); hold on;
plot(time, voltage2, 'b', 'LineWidth', 1.5);
% % xlabel('Time (s)');
% % ylabel('Voltage (V)');
% % title('Time vs Voltage');
% % legend('Voltage 1', 'Voltage 2');
% % xline(251)
xlim([1000, 2500])
ylim([-90  60])
% % box off
% % % % grid on;
% % 
% % % Define inset range (customize this!)
% % inset_start = stemTime - 15;   % e.g., start time in seconds
% % inset_end = stemTime + 25;   % e.g., end time in seconds
% % 
% % % Get indices for inset
% % inset_idx = time >= inset_start & time <= inset_end;
% % 
% % % Inset axes
% % inset_pos = [0.55 0.70 0.2 0.25];  % [x y width height] normalized
% % inset_axes = axes('Position', inset_pos);
% % box on;
% % plot(inset_axes, time(inset_idx), voltage1(inset_idx), 'r', 'LineWidth', 1.2); hold on;
% % plot(inset_axes, time(inset_idx), voltage2(inset_idx), 'b', 'LineWidth', 1.2);
% % % % title(inset_axes, sprintf('Zoom: %.1f–%.1f s', inset_start, inset_end));\
% % xlim([inset_start inset_end])
% % % % grid on;


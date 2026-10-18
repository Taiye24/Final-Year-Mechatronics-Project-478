%% generate_experimental_graphs.m
% Reproduces Figure 4.11 of Chapter 4: measured (a) drag force and
% (b) running trim against carriage speed for LCG 640, 560 and 480 mm.
% 640 mm: mean of three sets, error bars = u_x at 95% (random term, plus
% potentiometer resolution for trim;
% t95 = 4.303, N = 3). 560 and 480 mm: single sets.
% Data: ANALYZE DATA sheets of Data_Eval_640/560/480.xls (speeds are the measured run means).

clear; clc; close all;
outDir = fileparts(mfilename('fullpath'));

v640 = [0.41 0.92 1.44 1.96 2.47 2.98 3.50 4.02 4.54 5.05 5.56 6.08 6.59];
d640 = [1.06 4.97 16.48 28.59 32.26 34.47 37.16 40.79 43.28 45.40 47.29 49.71 52.49];
ud   = [0.56 0.18 0.84 0.07 0.11 1.10 0.40 0.48 0.32 0.70 0.64 0.23 0.77];
t640 = [2.25 2.32 2.67 3.93 4.69 4.95 5.01 5.33 5.58 5.38 5.08 4.75 4.45];
ut   = [0.15 0.18 0.27 0.05 0.10 0.18 0.13 0.12 0.08 0.08 0.09 0.08 0.09];  % includes b = 0.005 deg (half potentiometer resolution)

v560 = [0.39 0.90 1.42 1.94 2.45 2.97 3.48 4.00 4.51 5.02 5.54 6.06 6.57];
d560 = [0.77 5.43 15.96 31.12 37.85 41.39 45.95 47.71 47.79 48.21 48.70 49.99 51.38];
t560 = [3.25 3.43 3.96 5.70 6.45 6.78 7.26 7.47 7.06 6.50 5.95 5.47 5.02];

v480 = [0.39 0.90 1.42 1.93 2.45 2.92 3.48 4.00 4.51 5.03 5.54 6.05 6.57];
d480 = [1.30 5.91 16.80 41.24 48.63 52.26 57.52 55.96 53.43 51.56 50.87 51.77 50.59];
t480 = [4.40 4.61 5.44 7.69 8.53 9.45 9.88 9.16 8.29 7.51 6.76 6.42 5.68];

c640 = [0.122 0.306 0.475]; c560 = [0.773 0.353 0.067]; c480 = [0.329 0.510 0.208];
fig = figure('Color','w','Position',[100 100 1150 430]);

subplot(1,2,1); hold on;
h = errorbar(v640, d640, ud, '-o'); set(h, 'Color', c640, 'LineWidth', 1.6, 'MarkerSize', 5);
plot(v560, d560, '-s', 'Color', c560, 'LineWidth', 1.6, 'MarkerSize', 5);
plot(v480, d480, '-^', 'Color', c480, 'LineWidth', 1.6, 'MarkerSize', 5);
hold off; grid on;
xlabel('Speed, V (m/s)'); ylabel('Drag force, F (N)'); title('(a)');
legend('640 mm (0.40 L_p), 3 sets \pm u_x', ...
       '560 mm (0.35 L_p), 1 set', '480 mm (0.30 L_p), 1 set', ...
       'Location', 'southeast');

subplot(1,2,2); hold on;
h = errorbar(v640, t640, ut, '-o'); set(h, 'Color', c640, 'LineWidth', 1.6, 'MarkerSize', 5);
plot(v560, t560, '-s', 'Color', c560, 'LineWidth', 1.6, 'MarkerSize', 5);
plot(v480, t480, '-^', 'Color', c480, 'LineWidth', 1.6, 'MarkerSize', 5);
hold off; grid on;
xlabel('Speed, V (m/s)'); ylabel('Running trim, \tau (deg)'); title('(b)');

saveas(fig, fullfile(outDir, 'results_combined.png'));
fprintf('Saved results_combined.png to %s\n', outDir);

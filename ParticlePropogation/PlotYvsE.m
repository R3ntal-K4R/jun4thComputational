% Clear workspace, command window, and close all figures
clear; clc; close all;

% --- Data Section ---
% Define the path to the data file.
% This assumes there is a folder named 'output' in the same directory
% as this script, and inside 'output', there is a file '3_output.dat'.
% --- IMPORTANT: You can change this path to match your file ---
data_file_path = 'output/400_output.dat';

% Check if the file exists before trying to load it.
if exist(data_file_path, 'file')
    disp(['Loading data from: ' data_file_path]);
    % Use dlmread to load the data, skipping the first row (the header).
    % The '1' tells dlmread to start reading from the second row (row index 1).
    % The '02 tells it to start from the first column (column index 0).
    data = dlmread(data_file_path, '', 1, 0);
else
    % Display an error message and stop the script if the file doesn't exist.
    error('Data file not found at: %s\nPlease make sure the file exists and the path is correct.', data_file_path);
end

% --- Data Extraction ---
% In Octave, we can reference columns of a matrix easily.
% The ':' means "all rows".
Y = data(:, 4);  % The 'time' data is in the 2nd column
Eloss    = data(:, 13);  % The 'x' position data is in the 3rd column

% --- Plotting Section ---
% Create a new figure window.
figure;

% Plot the data. 'o-' creates a line plot with circular markers at each data point.
% 'LineWidth' and 'MarkerSize' are adjusted for better visibility.
plot(Y, Eloss, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0, 0.4, 0.7]);

% --- Formatting the Plot ---
% Add a title and labels to the axes.
title('Particle Trajectory: Y vs. Eloss');
xlabel('Y');
ylabel('Eloss');

% Add a grid for easier reading.
grid on;

% Improve the visual appearance of the plot
set(gca, 'FontSize', 12); % gca refers to the 'get current axes' handle

% Add a legend to identify the data series.
legend('Eloss', 'Location', 'northeast');

% Optional: Print a message to the command window
disp('Plot generated successfully.');

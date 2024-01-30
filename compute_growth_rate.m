% Author: Wolf De Wulf
% (but mostly co-pilot, this can probably be made much shorter)
%
% Example usage:
%
% [doubling_rate, doubling_time, growth_rate] = compute_growth_rate("data.txt", "OD_B", 5, "optimize", "optimize");
% --> this will optimise the smoothing method and number of slope change points to find the one that yields the highest R2
%
% BUT you can also manually specify them, atm only "gaussian" and "avg" are supported for smoothing:
% [doubling_rate, doubling_time, growth_rate] = compute_growth_rate("data.txt", "OD_B", 5, "gaussian", "optimize");
% [doubling_rate, doubling_time, growth_rate] = compute_growth_rate("data.txt", "OD_B", 5, "optimize", 2);
% --> note that you can choose the optimise one and specify the other
%
% TODO: - implement optmization of smoothing window
% TODO: - export values
%
% Changes by JB:
% - Fixed textbox behaviour
% - In Matlab the log function is the natural logarithm (LN).
% Therefore, linear regression to the LN(OD) curve will give us growth
% rate, from which we can calculate doubling time and then doubling rate.
% I amended the code accordingly. Alternatively, for the fit to give
% doubling rate directly, we would transform OD as follows: log10(smooth_od)/log10(2)

[doubling_rate, doubling_time, growth_rate] = compute_growth_ratee("data.txt", "OD_B", 5, "optimize", "optimize");
function [doubling_rate, doubling_time, growth_rate] = compute_growth_ratee(path, OD_var, smoothing_window_length, smoothing, max_num_changes)

% Read data
[od,time] = read_data(path, OD_var);

% Remove noisy values below reasonable detection threshold
od = od(od > 0.05);
time = time(od > 0.05);

plot_data(time, od, OD_var);

% Define search spaces
if string(smoothing) == "optimize"
    smoothing_methods = ["gaussian", "avg"];
else
    smoothing_methods = [smoothing];
end

if string(max_num_changes) == "optimize"
    max_num_changes_values = 2:10;
else
    max_num_changes_values = [max_num_changes];
end

% Search for the best smoothing method and max_num_changes
best_r2 = -Inf;
best_smoothing = "";
best_max_num_changes = 0;
for i = 1:length(smoothing_methods)
    for j = 1:length(max_num_changes_values)
        [~, r2, ~, ~, ~, ~, ~, ~] = compute_growth_rate_single(od, time, smoothing_methods(i), smoothing_window_length, max_num_changes_values(j));
        if r2 > best_r2
            best_r2 = r2;
            best_smoothing = smoothing_methods(i);
            best_max_num_changes = max_num_changes_values(j);
        end
    end
end

[smooth_od, r2, p, brkpt, od_fit, linear_time, linear_od, TF] = compute_growth_rate_single(od,time, best_smoothing, smoothing_window_length, best_max_num_changes);

% Compute growth rate, doubling rate, and doubling time
growth_rate = p(1);
doubling_time = log(2)*60/ growth_rate;
doubling_rate = (60/doubling_time);

% Plot
plot_growth_rate(OD_var, brkpt, linear_time, linear_od, od_fit, r2, doubling_rate, doubling_time, growth_rate, time, smooth_od, TF, p, best_max_num_changes, best_smoothing);

end

function [smooth_od, r2, p, brkpt, od_fit, linear_time, linear_od, TF] = compute_growth_rate_single(od, time, smoothing, smoothing_window_length, max_num_changes)
if smoothing == "gaussian"
    sigma = smoothing_window_length / 2; % Standard deviation for Gaussian filter
    filter_size = smoothing_window_length; % Size of the filter
    x = linspace(-filter_size / 2, filter_size / 2, filter_size);
    gaussian_filter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussian_filter = gaussian_filter' / sum(gaussian_filter); % Normalize
    smooth_od = conv(od, gaussian_filter, 'same');
elseif smoothing == "avg"
    smooth_od = movmean(od, smoothing_window_length);
end

% Detect the linear segment in log(y)
TF = ischange(log(smooth_od), 'linear', 'MaxNumChanges', max_num_changes);
brkpt = time(TF==1);

% Fit a linear regression to the linear segment
linear_time = time(time>=brkpt(1) & time<=brkpt(2));
linear_od = log(smooth_od(time>=brkpt(1) & time<=brkpt(2)));
p = polyfit(linear_time, linear_od, 1);

% Compute R^2
od_fit = polyval(p, linear_time);
residuals = linear_od - od_fit;
ss_res = sum(residuals.^2);
ss_tot = sum((linear_od - mean(linear_od)).^2);
r2 = 1 - ss_res / ss_tot;
end

function plot_data(time, smooth_od, OD_var)
figure; % Create a new figure window

% Plot the log od
plot(time, log(smooth_od));
title(sprintf('Growth Rate (%s)', OD_var));
xlabel('Time (hours)'); % Adjusted to reflect the change in time units
ylabel('LN(OD)');
end

function plot_growth_rate(OD_var, brkpt, linear_time, linear_od, od_fit, r2, doubling_rate, doubling_time, growth_rate, time, smooth_od, TF, p, max_num_changes, smoothing)
figure; % Create a new figure window

% Plot the log od
plot(time, log(smooth_od));
title(sprintf('Growth Rate (%s)', OD_var));
xlabel('Time (hours)'); % Adjusted to reflect the change in time units
ylabel('LN(OD)');
hold on; % Keep the current plot when adding new plots

% Plot the change points
plot(brkpt, log(smooth_od(TF==1)), 'ro');

% Plot the linear segment
plot(linear_time, linear_od, 'r', 'LineWidth', 2);

% Plot the linear regression line
od_fit = polyval(p, linear_time);
plot(linear_time, od_fit, 'k', 'LineWidth', 2);

% Display the equation, R-squared value, doubling rate, doubling time, and growth rate in the figure
str = sprintf('y = %.3fx + %.3f\nR^2 = %.3f\nDoubling Rate = %.3f 1/h\nDoubling Time = %.3f min\nGrowth Rate = %.3f 1/h\nMax num changes = %d\nSmoothing: %s', p(1), p(2), r2, doubling_rate, doubling_time, growth_rate, max_num_changes, smoothing);
% annotation('textbox', [0.78, 0.07, 0.2, 0.2], 'String', str, 'FitBoxToText', 'on' ,  'verticalalignment', 'bottom',  'horizontalalignment', 'right');
text(0.975,0.025, str, ...
    'units', 'normalized', ...
    'horizontalalignment','right', ...
    'verticalalignment', 'bottom');

% Add a legend
legend('log(y)','Slope Change Points', 'Linear Segment', 'Linear Fit');

% Set x-axis to be in hours and have ticks per hour
xlim([min(time), max(time)]);
xticks(min(time):1:max(time)); % Adjusted to reflect the change in time units
end

function [od, time] = read_data(path, OD_var)
% Open the file
fileID = fopen(path, 'r');

% Read the first line and split it into parts
firstLine = fgetl(fileID);
columns = split(firstLine, ',');
columns = columns(1:end-2); % Remove the last two entries

% Close the file
fclose(fileID);

% Read the data into a table
opts = detectImportOptions(path);
opts.VariableNames = columns;
opts.VariableTypes(1:numel(columns)) = {'double'}; % assuming all columns are numeric
data = readtable(path, opts);
od = data.(OD_var);
time = data.("time");
end
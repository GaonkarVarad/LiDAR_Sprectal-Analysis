%% MAIN SCRIPT: Spectra Shift Analysis
% Complete workflow for systematic error detection in spectral data

%% SECTION 1: DATA LOADING & VALIDATION
% Full absolute paths for MATLAB Online
dataA = readmatrix('/MATLAB Drive/Honeywell_Task/1d_A_0_5_tilt_EN.txt');
dataB = readmatrix('/MATLAB Drive/Honeywell_Task/1d_B_0_5_tilt_EN.txt');

% Extract coordinates with sanity check
if size(dataA,2) < 2 || size(dataB,2) < 2
    error('Data files must have at least 2 columns (X and Y values)');
end

xA = dataA(:,1); yA = dataA(:,2);
xB = dataB(:,1); yB = dataB(:,2);

%% SECTION 2: X² TRANSFORMATION WITH DUPLICATE HANDLING
[xA_sq, ia] = unique(xA.^2, 'stable'); 
yA_unique = yA(ia);  % Keep only unique X² values

[xB_sq, ib] = unique(xB.^2, 'stable');
yB_unique = yB(ib);

% Diagnostic output
fprintf('Unique points: A=%d/%d, B=%d/%d\n',...
        length(xA_sq), length(xA), length(xB_sq), length(xB));

%% SECTION 3: COMMON GRID & INTERPOLATION
grid_min = max(min(xA_sq), min(xB_sq));
grid_max = min(max(xA_sq), max(xB_sq));
common_grid = linspace(grid_min, grid_max, 1000)';

% Safe interpolation (pchip handles monotonic data well)
yA_interp = interp1(xA_sq, yA_unique, common_grid, 'pchip');
yB_interp = interp1(xB_sq, yB_unique, common_grid, 'pchip');

%% SECTION 4: SHIFT CALCULATION (CROSS-CORRELATION)
[corr_vals, lags] = xcorr(yA_interp, yB_interp);
[~, max_idx] = max(corr_vals);
optimal_lag = lags(max_idx);

% Convert lag to X units
dx_sq = mean(diff(common_grid));
x_shift = sqrt(common_grid(1) + optimal_lag*dx_sq) - sqrt(common_grid(1));

%% SECTION 5: VALIDATION & VISUALIZATION
figure;

% Subplot 1: Original Data Comparison
subplot(3,1,1);
plot(xA, yA, 'b', xB, yB, 'r');
title('Raw Spectra Comparison');
legend('Spectrum A', 'Spectrum B');

% Subplot 2: X² Transformed Data
subplot(3,1,2);
plot(common_grid, yA_interp, 'b', common_grid, yB_interp, 'r');
title('X² Transformed Spectra');
legend('Transformed A', 'Transformed B');

% Subplot 3: Alignment Result
subplot(3,1,3);
plot(common_grid, yA_interp, 'b',...
     common_grid + optimal_lag*dx_sq, yB_interp, 'r--');
title('Aligned Spectra');
legend('Spectrum A', 'Shifted Spectrum B');

%% SECTION 6: RESULTS OUTPUT
fprintf('\n=== ANALYSIS RESULTS ===\n');
fprintf('Optimal Lag: %d samples\n', optimal_lag);
fprintf('Systematic X-shift: %.4f units\n', x_shift);
fprintf('Grid Resolution: %.4e units/sample\n', dx_sq);

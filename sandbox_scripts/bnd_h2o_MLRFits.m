function [fitresult, gof] = bnd_h2o_MLRFits %(Q,pf)
%% MULTILINEAR FIT FUNCTIONS for BND and H20_ex
% These are crude esimators only.
% (BND is the more important)
%  Input:  Q  = mass discharge rate (kg/s)
%          pf = ambient pressure at vent (Pa)
%
%  Output: BND = bubble number density (m^-3) after first nucleation
%          h20_ex = mass fraction of exsolved h20 after first nucleation
%                   (ie h20_ex/h20_0)
%
% CREATEFITS(X,Y,Z,Z2)
%  Create fits.
%
%  Data for 'BND = f(Q,pf)' fit:
%      X Input : x
%      Y Input : y
%      Z Output: z
%  Data for 'H2O_ex = f(Q,pf)' fit:
%      X Input : x
%      Y Input : y
%      Z Output: z2
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Dec-2023 13:41:14

fitFile = 'BND_&_H20_ex_multilinear_fits';
load(fitFile)
% 
% Q = log10(Q);
%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 2, 1 );
gof = struct( 'sse', cell( 2, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'BND = f(Q,pf)'.
% CRUDE predictor for Bubble Number Density after first nucleation
%  `-> function of MER and surface pressure
%  `-> uses numerical results of Rowell et al 2022 to make this prediction
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = fittype( 'poly11' );

% Fit model to data.
[fitresult{1}, gof(1)] = fit( [xData, yData], zData, ft );

% Plot fit with data.
figure( 'Name', 'BND = f(Q,pf)' );
h = plot( fitresult{1}, [xData, yData], zData );
legend( h, 'BND = f(Q,pf)', 'z vs. x, y', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
zlabel( 'z', 'Interpreter', 'none' );
grid on
view( -12.3, 9.8 );

%% Fit: 'H2O_ex = f(Q,pf)'.
% CRUDE predictor for exsolved fraction of h20 after first nucleation
%  `-> function of MER and surface pressure
%  `-> uses numerical results of Rowell et al 2022 to make this prediction
[xData, yData, zData] = prepareSurfaceData( x, y, z2 );

% Set up fittype and options.
ft = fittype( 'poly11' );

% Fit model to data.
[fitresult{2}, gof(2)] = fit( [xData, yData], zData, ft );

% Plot fit with data.
figure( 'Name', 'H2O_ex = f(Q,pf)' );
h = plot( fitresult{2}, [xData, yData], zData );
legend( h, 'H2O_ex = f(Q,pf)', 'z2 vs. x, y', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
zlabel( 'z2', 'Interpreter', 'none' );
grid on
view( -44.5, 17.7 );
end
% MATLAB Code to Plot Electron Cloud Weight Map for p_z Orbital with Sign-based Coloring

% Clear workspace and command window
clear;
clc;

% Define parameters
a0 = 0.01; % Bohr radius (arbitrary units)
rMax = 2; % Maximum radius
numPoints = 50; % Number of points in each direction

% Create spherical grid
r = linspace(0, rMax, numPoints);
theta = linspace(0, pi, numPoints);
phi = linspace(0, 2*pi, numPoints);

% Create meshgrid for spherical coordinates
[R, Theta, Phi] = ndgrid(r, theta, phi);

% Define the radial function R(r)
R_r = sqrt(2/(a0^3)) .* (R/a0) .* exp(-R/(2*a0));

% Define the angular function Y(theta)
Y_theta = sqrt(5/(4*pi)) .* cos(Theta);

% Calculate the p_z orbital wave function
psi_pz = R_r .* Y_theta;

% Calculate probability density
probabilityDensity = abs(psi_pz).^2;

% Convert spherical coordinates to Cartesian coordinates for plotting
X = R .* sin(Theta) .* cos(Phi);
Y = R .* sin(Theta) .* sin(Phi);
Z = R .* cos(Theta);

% Create a 3D figure for the weight map
figure;

% Define threshold value for isosurface
threshold = 0.01; % Use a low threshold to show all data

% Plot positive isosurface
idx1 = psi_pz<0;
pb1 = probabilityDensity;
pb1(idx1) = 0;
p_pos = patch(isosurface(X, Y, Z, pb1, threshold));
set(p_pos, 'FaceColor', 'magenta', 'EdgeColor', 'none', 'FaceAlpha', 0.7); 
hold on;

% Plot negative isosurface with a different color
% Create a negative probability density for negative values
idx2 = psi_pz>0;
pb2 = probabilityDensity;
pb2(idx2) = 0;
p_neg = patch(isosurface(X, Y, Z, pb2, threshold)); % Set the threshold to negative
set(p_neg, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.7); 

% Set view and labels
view(3);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Electron Cloud Weight Map of p_z Orbital with Sign-based Coloring');
axis equal;
grid on;

% Optional: Add lighting for better visualization
camlight;
lighting gouraud;

% Hold off to stop adding to current figure
hold off;

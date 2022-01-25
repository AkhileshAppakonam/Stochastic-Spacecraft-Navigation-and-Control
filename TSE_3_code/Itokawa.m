load('ItokawaData.mat')
load('ItokawaIndex.mat')
load('ItokawaIndex_plus1.mat')

% A is the xyz matrix for Pan data
% B is the matrix of Pan's indices indexed from 0
% C is the matrix of Pan's indices indexed from 1

x = I(1:length(I), 1);
y = I(1:length(I), 2);
z = I(1:length(I), 3);

% Creates a surface and applies texture (increase matrix
% dimesnsions for N & M to increase resolution of graph)
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

N = 5000;  % Number of y values in uniform grid
M = 5000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x(:), y(:), z(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)
image = imread('Itokawa.bmp');             % Load RGB image

figure(1)
h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');
grid on;
axis equal;
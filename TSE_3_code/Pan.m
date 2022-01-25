load('PanData.mat')
load('PanIndex.mat')
load('PanIndex_plus1.mat')

% A is the xyz matrix for Pan data
% B is the matrix of Pan's indices indexed from 0
% C is the matrix of Pan's indices indexed from 1

x = A(1:length(A), 1);
y = A(1:length(A), 2);
z = A(1:length(A), 3);


% Creates a surface and applies texture (increase matrix
% dimesnsions for N & M to increase resolution of graph)
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

N = 2000;  % Number of y values in uniform grid
M = 2000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x(:), y(:), z(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)
image = imread('Pan.bmp');             % Load RGB image

figure(1)
h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');
grid on;
axis equal;

%% 

% This section of code is to obtain a high resolution graph using minimal matrice dimensions by splitting the model into its top and bottom halves.
% Run this section if you are unable to run the code above with values of N and M equal to or above 10000 points

load('PanData.mat')
load('PanIndex.mat')
load('PanIndex_plus1.mat')

x = A(1:length(A), 1);
y = A(1:length(A), 2);
z = A(1:length(A), 3);

x_bot = A(1:6633, 1);
y_bot = A(1:6633, 2);
z_bot = A(1:6633, 3);

x_top = A(6634:length(A), 1);
y_top = A(6634:length(A), 2);
z_top = A(6634:length(A), 3);

image = imread('Pan.bmp');

% Creates bottom surface and applies texture 
% High res. can be obtained with N-by-M dimensions of just 1000-by-1000
xmin = min(x_bot);
xmax = max(x_bot);
ymin = min(y_bot);
ymax = max(y_bot);

N = 3000;  % Number of y values in uniform grid
M = 3000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x_bot(:), y_bot(:), z_bot(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)


figure(2)
h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');
grid on;
axis equal;
hold on;


% Creates top surface and applies texture 
% High res. can be obtained with N-by-M dimensions of just 1000-by-1000
xmin = min(x_top);
xmax = max(x_top);
ymin = min(y_top);
ymax = max(y_top);

N = 3000;  % Number of y values in uniform grid
M = 3000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x_top(:), y_top(:), z_top(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)

h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');

     
% Fills in missing middle portion of the model which connects the two halves     
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

N = 3000;  % Number of y values in uniform grid
M = 3000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x(:), y(:), z(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)
image = imread('Pan.bmp');             % Load RGB image

h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');   

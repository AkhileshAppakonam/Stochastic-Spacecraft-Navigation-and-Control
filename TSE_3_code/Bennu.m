load('BennuData.mat')
load('BennuIndex.mat')
load('BennuIndex_plus1.mat')

% A is the xyz matrix for Bennu data
% B is the matrix of Bennu indices indexed from 0
% C is the matrix of Bennu indices indexed from 1

x = A(1:1148, 1);
y = A(1:1148, 2);
z = A(1:1148, 3);

% Creates a surface and applies texture (increase matrix
% dimesnsions for N & M to increase resolution of graph)
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
image = imread('Bennu.bmp');             % Load RGB image

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

load('BennuData.mat')
load('BennuIndex.mat')
load('BennuIndex_plus1.mat')

x = A(1:1148, 1);
y = A(1:1148, 2);
z = A(1:1148, 3);

x_top = A(1:574, 1);
y_top = A(1:574, 2);
z_top = A(1:574, 3);

x_bot = A(575:1148, 1);
y_bot = A(575:1148, 2);
z_bot = A(575:1148, 3);

image = imread('Bennu.bmp');

% Creates bottom surface and applies texture 
% High res. can be obtained with N-by-M dimensions of just 1000-by-1000
xmin = min(x_bot);
xmax = max(x_bot);
ymin = min(y_bot);
ymax = max(y_bot);

N = 1000;  % Number of y values in uniform grid
M = 1000;  % Number of x values in uniform grid
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

N = 1000;  % Number of y values in uniform grid
M = 1000;  % Number of x values in uniform grid
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

N = 2000;  % Number of y values in uniform grid
M = 2000;  % Number of x values in uniform grid
xu = linspace(xmin, xmax, M);         % Uniform x-coordinates
yu = linspace(ymin, ymax, N);         % Uniform y-coordinates
[X, Y] = meshgrid(xu, yu);            % Create meshes for xu and yu
F = TriScatteredInterp(x(:), y(:), z(:));  % Create interpolant
Z = F(X, Y);                          % Evaluate interpolant (N-by-M matrix)
image = imread('Bennu.bmp');             % Load RGB image

h = surf(X, Y, Z);
set(h, 'Cdata', flip(image, 1), ...  % Plot surface
         'FaceColor', 'texturemap', ...
         'EdgeColor', 'none');   

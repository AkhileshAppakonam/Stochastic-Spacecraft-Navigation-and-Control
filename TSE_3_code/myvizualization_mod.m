function myvizualization_mod(t,B,A )

global cont_choose pw tumbling Asteroid

if Asteroid== "Pan"
Lx = 6;     Ly = 9;     Lz = 3;
R_Asteroid = 30;
World_size = 100;

elseif Asteroid== "Bennu"
Lx = 6/100;     Ly = 9/100;     Lz = 3/100;
R_Asteroid = 0.536;
World_size = 2;
end

% Block size (length of the block along the x,y,z coodinates)
% Lx = 6;     Ly = 9;     Lz = 3;
%Lx = 2;     Ly = 3;     Lz = 1;


% The vizualization will be +/- World_size in all
% three axes directions (x,y,z)

% % ----------- Motion data  -----------
% t :  Time data size [t by 1]
%
% r :  Position data w.r.t. the origin of the vizualization volume
%   :  The position data of the four blocks is stored
%   :  in a matrix size [t by 12], where t is all the Time steps data
%   :  r=[ r1 r2 r3 r4]=[ 1:3 4:6 7:9 10:12]
%
%
%
% A :  Orientation data (x-y-z Euler angle) about block C.M.
%   :  The rientation data of the four blocks is stored
%   :  in a matrix size [t by 12], where t is all the Time steps data
%   :  A=[ A1 A2 A3 A4]=[ 1:3 4:6 7:9 10:12]

% Extract the total number of time steps in order to generate
% the vizualization
n_time = length(t);

% Compute propagation of vertices and patches
for i_time=1:n_time
    %     R = EulerR(A(i_time,1:3));
    R=A(1:3,1:3,i_time);
    r=B(:,:,i_time);
    VertexData(:,:,i_time) = VizVerMakeBlock(r,R,[Lx,Ly,Lz]);
    [X,Y,Z] = VizPatMakeBlock(VertexData(:,:,i_time));
    PatchData_X(:,:,i_time) = X;
    PatchData_Y(:,:,i_time) = Y;
    PatchData_Z(:,:,i_time) = Z;
end

% Draw initial figure
h = patch(PatchData_X(:,:,1),PatchData_Y(:,:,1),PatchData_Z(:,:,1),'y');
set(h,'FaceLighting','phong','EdgeLighting','phong');

hold on, plot3(squeeze(B(1,1,:)),squeeze(B(1,2,:)),squeeze(B(1,3,:)),'k--')
[X,Y,Z] = sphere(25);
surface(X*R_Asteroid,Y*R_Asteroid,Z*R_Asteroid)

% Draw Title, Labels, and Inital camera angle
%title(['Time = ', num2str(t(j))]);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
set(gca,'FontSize',14);
axis vis3d equal;
view([-37,30]);
camlight;
grid on;

% Stablish the size of vizualization volume
xlim([-World_size,World_size]);
ylim([-World_size,World_size]);
zlim([-World_size,World_size]);

% -----------------------  Set up the movie. ----------------------
% -----------------------------------------------------------------
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 30;           % How many frames per second.
% -----------------------------------------------------------------

open(writerObj);

for i_time=1:n_time
    
    set(h,'XData',PatchData_X(:,:,i_time));
    set(h,'YData',PatchData_Y(:,:,i_time));
    set(h,'ZData',PatchData_Z(:,:,i_time));
    
    if strcmp(cont_choose, 'present work') == 1 && ...
            strcmp(pw, 'TSE(3)') == 1 && strcmp(tumbling, 'tumbling') == 1
        pause(.02);
    elseif strcmp(tumbling, 'no tumbling') == 1
        pause(.015);
    else
        %pause(.0005);
        % pause(.5);
    end
    drawnow;
    frame = getframe(gcf);        % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame); % Saves the movie.    
end
close(writerObj);                 % Close the movie file.

end


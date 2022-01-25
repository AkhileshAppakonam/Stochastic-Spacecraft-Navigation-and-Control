%===================================================================================================
%[]FUNCTION NAME: AnimateGroundTrack.m
%[]AUTHOR: Julio César Benavides, Revised by Matthew Wittal
%[]CREATED: 06/24/2010
%[]REVISED: 07/23/2020
%===================================================================================================
%[]FUNCTION DESCRIPTION:
%This function animates a satellite ground track with respect to the the
%given asteroid.
%===================================================================================================
%[]INPUT VARIABLE:
%(t)|Simulation time:
%---------------------------------------------------------------------------------------------------
%(R)|Satelllite position vectors with respect to the Earth in Earth-Centered Inertial coordinates.
%---------------------------------------------------------------------------------------------------
%(UTC)|Coordinated universal time at mission commencement.
%---------------------------------------------------------------------------------------------------
%(dt)|Time step.
%===================================================================================================
%[]OUTPUT VARIABLES:
%None.
%===================================================================================================
%[]VARIABLE FORMAT:
%(t)|Column Vector {m x 1} or Row Vector {1 x n}.
%---------------------------------------------------------------------------------------------------
%(R)|Matrix {3 x n}.
%---------------------------------------------------------------------------------------------------
%(UTC)|Row Vector {1 x 6}.
%---------------------------------------------------------------------------------------------------
%(dt)|Scalar {1 x 1}.
%===================================================================================================
%[]AUXILIARY FUNCTIONS:
%(ROTATION.m)|This function returns the rotation matrix for a given angle about a given axis.
%===================================================================================================
%[]COMMENTS:
%This function requires the image file "GROUNDTRACK.bmp" in order to function properly.
%===================================================================================================
function AnimateGroundTrack_New(t,R,R2,dt,Asteroid)

    if Asteroid == "Pan"
        
        Wp = 2 * pi / 49684.3812 * [0; 0; 1];
        %[rad/s]Angular velocity of Pan in P.C.E.F. coordinates.
        
        Implot = imread('Pan.bmp');
        %[]Loads the groundtrack image data.
        
    elseif Asteroid == "Bennu"
        
        Wp = 2 * pi / (4.297*60*60) * [0; 0; 1];
        %[rad/s]Angular velocity of Bennu in P.C.E.F. coordinates.
        
        Implot = imread('Bennu.bmp');
        %[]Loads the groundtrack image data.
        
    elseif Asteroid == "Itokawa"
        
        Wp  = 2 * pi / (12.132*60*60) * [0; 0; 1];
        %[rad/s]Angular velocity of Itokawa in P.C.E.F. coordinates.
        
        Implot = imread('Itokawa.bmp');
        %[]Loads the groundtrack image data.
        
    end
    
    PRA = wrapTo2Pi(norm(Wp) * t);
    %[rad]Earth rotation angle at mission commencement.
    
    s = numel(t);
    %[]Number of satellite positions that are being plotted.
    
    latitude = zeros(1,s);
    %[]Allocates memory for the satellite latitudes during the simulation time.
    
    longitude = zeros(1,s);
    %[]Allocates memory for the satellite longitudes during the simulation time.
    
    latitude2 = zeros(1,s);
    %[]Allocates memory for the satellite latitudes during the simulation time.
    
    longitude2 = zeros(1,s);
    %[]Allocates memory for the satellite longitudes during the simulation time.
    
    for k = 1:s
        
        EciToEcef = [cos(PRA(k)), sin(PRA(k)), 0; -sin(PRA(k)), cos(PRA(k)), 0; 0, 0, 1];
        %[]Current matrix that transforms vectors from ECI coordinates to ECEF coordinates.
        
        Recef = EciToEcef * R(1:3,k);
        %[km]Current satellite position WRT the Earth in ECEF coordinates.
        
        Recef2 = EciToEcef * R2(1:3,k);
        %[km]Current satellite position WRT the Earth in ECEF coordinates.
        
        latitude(k) = asin(Recef(3) / norm(Recef)) * 180 / pi;
        %[deg]Current satellite latitude.
        
        longitude(k) = atan2(Recef(2),Recef(1)) * 180 / pi;
        %[deg]Current satellite longitude.
        
        latitude2(k) = asin(Recef2(3) / norm(Recef2)) * 180 / pi;
        %[deg]Current satellite latitude.
        
        longitude2(k) = atan2(Recef2(2),Recef2(1)) * 180 / pi;
        %[deg]Current satellite longitude.
        
    end
    
    Window2 = figure(3);
    %[]Opens a new window.
    
    GroundTrack = image(Implot);
    %[]Adds a plot to the current axes.
    
    Axes = axes('Parent',Window2);
    %[]Adds an axes to the specified window.
    
    ScreenSize = get(0,'ScreenSize');
    %[]Determines the location and dimensions of the current monitor.
    
    set(GroundTrack, ...
        'XData',[-180,180], ...
        'YData',[-90,90], ...
        'Parent',Axes);
    %[]Adjusts the properties of the groundtrack image.
    
    set(Window2, ...
        'NumberTitle','Off', ...
        'Name','GROUND TRACK', ...
        'Color','k',...
        'OuterPosition',ScreenSize);
    %[]Adjusts the properties of the specified window.
    
    string =    {'180W' [-90]
                '165W'  [-82.5]
                '150W'  [-75]
                '135W'  [-67.5]
                '120W'  [-60]
                '105W'  [-52.5]
                '90W'   [-45]
                '75W'   [-37.5]
                '60W'   [-30]
                '45W'   [-22.5]
                '30W'   [-15]
                '15W'   [-7.5]
                '0'     [0]
                '15E'   [7.5]
                '30E'   [15]
                '45E'   [22.5]
                '60E'   [30]
                '75E'   [37.5]
                '90E'   [45]
                '105E'  [52.5]
                '120E'  [60]
                '135E'  [67.5]
                '150E'  [75]
                '165E'  [82.5]
                '180E'  [90]};
    %[]Longitude string.
    
    string2 =   {'90N'  [-180]
                '75N'   [-150]
                '60N'   [-120]
                '45N'   [-90]
                '30N'   [-60]
                '15N'   [-30]
                '0'     [0]
                '15S'   [30]
                '30S'   [60]
                '45S'   [90]
                '60S'   [120]
                '75S'   [150]
                '90S'   [180]};
    %[]Latitude string
    
    set(Axes, ...
        'FontName','Arial', ...
        'FontSize',8, ...
        'FontWeight','Bold', ...
        'NextPlot','Add', ...
        'XColor','w', ...
        'YColor','w', ...
        'XGrid','On', ...
        'YGrid','On', ...
        'XLim',[-180,180], ...
        'YLim',[-90,90], ...
        'XTick',-180:15:180, ...
        'YTick',-90:15:90, ...
        'XTickLabel',string, ...
        'YTickLabel',string2);
    %[]Adjusts the properties of the specified axis.
    
    xlabel('Longitude (\circ)','FontSize',12,'Parent',Axes);
    %[]Adds a label to the specified x-axis and adjusts its properties.
    
    ylabel('Latitude (\circ)','FontSize',12,'Parent',Axes);
    %[]Adds a label to the specified y-axis and adjusts its properties.
    
    counter = 0;
    %[s]Clock counter.
    
    while counter <= max(t)
        
        index = find(t <= counter,1,'Last');
        %[]Index for the new time.
        
        plot(longitude(index),-latitude(index), ...
            'Color','r', ...
            'LineStyle','None', ...
            'Marker','.', ...
            'MarkerSize',15, ...
            'Parent',Axes);
        %[]Plots the satellite ground track.
        
        plot(longitude2(index),-latitude2(index), ...
            'Color','b', ...
            'LineStyle','None', ...
            'Marker','.', ...
            'MarkerSize',15, ...
            'Parent',Axes);
        %[]Plots the satellite ground track.
        
        counter = counter + dt;
        %[s]Clock counter update.
        
        pause(0.1);
        %[s]Pauses the loop for the specified amount of time.
        
    end
    
    title('Simulation Complete!','Color','w','FontSize',12,'Parent',Axes);
    %[]Adds a title to the specified axes and adjusts its properties.
    
end
%===================================================================================================
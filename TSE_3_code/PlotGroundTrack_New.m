%===================================================================================================
%[]FUNCTION NAME: PlotGroundTrack.m
%[]AUTHOR: Julio Cesar Benavides, Revised by Matthew Wittal
%[]CREATED: 01/15/2013
%[]REVISED: 07/23/2020
%===================================================================================================
%[]FUNCTION DESCRIPTION:
%This function plots a satellite's ground track with respect to the a
%number of given asteroids.
%===================================================================================================
%[]INPUT VARIABLES:
%(t)|Time vector.
%---------------------------------------------------------------------------------------------------
%(R)|Satelllite position vectors with respect to the Earth in Earth-Centered Inertial coordinates.
%---------------------------------------------------------------------------------------------------
%(UTC)|Coordinated universal time at mission commencement.
%===================================================================================================
%[]OUTPUT VARIABLES:
%None.
%===================================================================================================
%[]VARIABLE FORMAT:
%(t)|Column Vector {m x 1} or Row Vector {1 x n}.
%---------------------------------------------------------------------------------------------------
%(R * R2)|Matrix {3 x n}.
%---------------------------------------------------------------------------------------------------
%(UTC)|Row Vector {1 x 6}.
%===================================================================================================
%[]USER-DEFINED FUNCTIONS:
%None.
%===================================================================================================
%[]COMMENTS:
%The universal time must be a row vector defined as UTC = [year, month, day, hour, minute, second].
%'WRT' stands for "with respect to".  'ECI' stands for "Earth-Centered Inertial".  'ECEF' stands for
%"Earth-Centered Earth-Fixed".
%===================================================================================================
function PlotGroundTrack_New(t,R,R2,Asteroid)
    
    if Asteroid == "Pan"
        
        Wp = 2 * pi / 49684.3812 * [0; 0; 1];
        %[rad/s]Angular velocity of Pan in P.C.E.F. coordinates.
        
        Implot = imread('Pan.bmp');
        %[]Loads the groundtrack image data.
        
        GroundTrack = image(Implot);
        %[]Adds a plot to the current axes.
        
    elseif Asteroid == "Bennu"
        
        Wp = 2 * pi / (4.297*60*60) * [0; 0; 1];
        %[rad/s]Angular velocity of Bennu in P.C.E.F. coordinates.
        
        Implot = imread('Bennu.bmp');
        %[]Loads the groundtrack image data.
        
        GroundTrack = image(Implot);
        %[]Adds a plot to the current axes.
        
    elseif Asteroid == "Itokawa"
        
        Wp  = 2 * pi / (12.132*60*60) * [0; 0; 1];
        %[rad/s]Angular velocity of Itokawa in P.C.E.F. coordinates.
        
        Implot = imread('Itokawa.bmp');
        %[]Loads the groundtrack image data.
        
        GroundTrack = image(Implot);
        %[]Adds a plot to the current axes.
        
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
    
    Window1 = figure(3);
    %[]Opens a new window.
    
    Axes = axes('Parent',Window1);
    %[]Adds an axes to the specified window.
    
    ScreenSize = get(0,'ScreenSize');
    %[]Determines the location and dimensions of the current monitor.
    
    set(GroundTrack, ...
        'XData',[-180,180], ...
        'YData',[-90,90], ...
        'Parent',Axes);
    %[]Adjusts the properties of the groundtrack image.
    
    set(Window1, ...
        'Color','w', ...
        'NumberTitle','Off', ...
        'Name','SATELLITE GROUND TRACK', ...
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
        'XColor','k', ...
        'YColor','k', ...
        'XGrid','On', ...
        'YGrid','On', ...
        'XLim',[-180,180], ...
        'YLim',[-90,90], ...
        'XTick',-180:15:180, ...
        'YTick',-90:15:90, ...
        'XTickLabel',string, ...
        'YTickLabel',string2);
    %[]Adjusts the properties of the specified axis.
    
    plot(longitude,-latitude, ...
        'Color','r', ...
        'LineStyle','None', ...
        'Marker','.', ...
        'MarkerSize',10, ...
        'Parent',Axes);
    %[]Plots the satellite ground track.
   
    
    plot(longitude2,-latitude2, ...
        'Color','b', ...
        'LineStyle','None', ...
        'Marker','.', ...
        'MarkerSize',10, ...
        'Parent',Axes);
    %[]Plots the satellite ground track.
    
    xlabel('Longitude (\circ)','FontSize',12,'Parent',Axes);
    %[]Adds a label to the x-axis.
    
    ylabel('Latitude (\circ)','FontSize',12,'Parent',Axes);
    %[]Adds a label to the y-axis.
    
end
%===================================================================================================
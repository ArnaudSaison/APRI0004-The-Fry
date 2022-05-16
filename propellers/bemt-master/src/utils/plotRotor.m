function plotRotor(Blades, Elem)
% PLOTROTOR Display the blade and rotor geometries, using the airfoil specified in CONFIG_FILE.
% -----
%
% Synopsis: [-] = PLOTROTOR(Blades, Elem)
%
% Input:    Blades = (required) Blade parameters defined in config file
%           Elem   = (required) Every data related to blade elements (chords, pitch,...)
%
% Output:   2 Figures (Single Blade and Full Rotor)
%
% Calls:    none
%
% See also: BEMT.

% ------
% This script is an optional part of the BEMT code.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------

%% Parameters
SECTION_SPACING = 1;    % Plot lines every N sections on the blades


%% Load and prepare coordinate systems

% Load airfoil coordinates
profileFileID = fopen(['./airfoils_data/',Blades.PROFILE_FILE]);
airfoilCoord = cell2mat(textscan(profileFileID,'%f %f'));
fclose(profileFileID);

% Manipulate coordinates to form a continuous loop if required
isLoop = airfoilCoord(1,1) == airfoilCoord(end,1);
if ~isLoop
    airfoilCoord(end/2+1:end,:)=flip(airfoilCoord(end/2+1:end,:));
end


% Chord of each element
chordMatrix = repmat(Elem.Chord,length(airfoilCoord(:,1)),1);
thrustMatrix = repmat(Elem.dT,length(airfoilCoord(:,1)),1);


% Scale and position of the airfoil sections
initGeomX=repmat(airfoilCoord(:,1)-0.5,1,Blades.nELEM).*chordMatrix ;
initGeomY=repmat(Elem.Ypos,length(airfoilCoord(:,1)),1);
initGeomZ=repmat(airfoilCoord(:,2),1,Blades.nELEM).*chordMatrix;

% Pitch the airfoil sections according to the local angle (collective + twist)
for i=1:Blades.nELEM
    RotY=[cos(Elem.Pitch(i)) 0 sin(Elem.Pitch(i));0 1 0; -sin(Elem.Pitch(i)) 0 cos(Elem.Pitch(i))];
    for j=1:length(airfoilCoord(:,1))
        dummy=RotY*[initGeomX(j,i); initGeomY(j,i); initGeomZ(j,i)];
        geomX(j,i)=dummy(1);
        geomY(j,i)=dummy(2);
        geomZ(j,i)=dummy(3);
    end
end


% Repeat the blade around the center to create the full rotor
RotBl(1).X=geomX;
RotBl(1).Y=geomY;
RotBl(1).Z=geomZ;

for iBl = 2:Blades.nBLADES
    bladePosAngle=2*pi/Blades.nBLADES*(iBl-1);
    RotZ=[cos(bladePosAngle) -sin(bladePosAngle) 0;sin(bladePosAngle) cos(bladePosAngle) 0; 0 0 1];
    for i=1:length(geomX(:,1))
        for j=1:length(geomX(1,:))
            dummy=RotZ*[RotBl(1).X(i,j); RotBl(1).Y(i,j); RotBl(1).Z(i,j)];
            tmpGeomX(i,j)=dummy(1);
            tmpGeomY(i,j)=dummy(2);
            tmpGeomZ(i,j)=dummy(3);
        end
    end
    RotBl(iBl).X=tmpGeomX;
    RotBl(iBl).Y=tmpGeomY;
    RotBl(iBl).Z=tmpGeomZ;
end

% Create a object for the rotor hub
[hubX, hubY, hubZ]=cylinder(Blades.ROOT_CUTOUT*1.25,18);
hubZ(1,:)=max(max(RotBl(1).Z))*1.25;
hubZ(2,:)=min(min(RotBl(1).Z))*1.25;

% Create a cone for top of hub
r = -linspace(0,Blades.ROOT_CUTOUT*1.25,2) ;
th = linspace(0,2*pi,19) ;
[R,T] = meshgrid(r,th) ;
coneX = R.*cos(T) ;
coneY = R.*sin(T) ;
coneZ = 0.5*R+hubZ(1,1)+0.5*(R(1,1)-R(1,end)) ;


%% Plot the blade and rotor

% Single Blade
surf2stl('blade.stl', geomX,geomY,geomZ)

figname='Single Blade';
figure('PaperUnits', 'inches', 'PaperPosition', [0 0 1280 1024]/250,'Name',figname)
set(gcf,'units','points','position',[0,1000,800,600])

surf(geomX,geomY,geomZ,thrustMatrix,'EdgeColor','none','Linestyle','-')
hold on
plot3(geomX(:,1:SECTION_SPACING:end),geomY(:,1:SECTION_SPACING:end),geomZ(:,1:SECTION_SPACING:end),'-k')
hold off
grid on
axis equal
hCol=colorbar;
hCol.Label.String = 'Sectional thrust, dT [N]';
set(gca,...
    'Box', 'off',...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XGrid'       , 'on'      , ...
    'XColor'      , 'k', ...
    'YColor'      , 'k', ...
    'GridLineStyle', ':'         , ...
    'GridColor', 'k', ...
    'GridAlpha', 0.25, ...
    'LineWidth'   , 1         );
set(gcf, 'PaperPositionMode', 'auto');

% Full Rotor
figname='Full Rotor';
figure('PaperUnits', 'inches', 'PaperPosition', [0 0 1280 1024]/250,'Name',figname)
set(gcf,'units','points','position',[800,1000,800,600])

surf(RotBl(1).X,RotBl(1).Y,RotBl(1).Z,thrustMatrix,'EdgeColor','none')
hold on
plot3(RotBl(1).X(:,1:SECTION_SPACING:end),RotBl(1).Y(:,1:SECTION_SPACING:end),RotBl(1).Z(:,1:SECTION_SPACING:end),'-k')
for iBl=2:Blades.nBLADES
    surf(RotBl(iBl).X,RotBl(iBl).Y,RotBl(iBl).Z,thrustMatrix,'EdgeColor','none')
    plot3(RotBl(iBl).X(:,1:SECTION_SPACING:end),RotBl(iBl).Y(:,1:SECTION_SPACING:end),RotBl(iBl).Z(:,1:SECTION_SPACING:end),'-k')
end
surf(hubX,hubY,hubZ,'FaceColor',[0 0.4470 0.7410])
surf(coneX,coneY,coneZ,'FaceColor',[0 0.4470 0.7410])
axis equal
grid on
hCol=colorbar;
hCol.Label.String = 'Sectional thrust, dT [N]';
set(gca,...
    'Box', 'off',...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'on'      , ...
    'XGrid'       , 'on'      , ...
    'XColor'      , 'k', ...
    'YColor'      , 'k', ...
    'GridLineStyle', ':'         , ...
    'GridColor', 'k', ...
    'GridAlpha', 0.25, ...
    'LineWidth'   , 1         );
set(gcf, 'PaperPositionMode', 'auto');


end

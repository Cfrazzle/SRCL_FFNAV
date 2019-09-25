%% =============================Animate================================= %%
%
%                             DESCRIPTION:
% This routine generates the animations for a given data set
%
close all
global rE ERR_r

%% Position Error 

%Initializing the Drawing Space
set(gcf,'Renderer','zbuffer','Menubar','default','Name','Satellite Formation Simulation', ... 
    'NumberTitle','off','Position',[300,50,650,650], ... 
    'Color',[01 01 01]); 
lim=ERR_r;%Setting the limits of the graph
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
title('Satellite Formation Simulation');


v = VideoWriter([pwd '\animations\Error_.mpeg']); 

open(v)
for i=1:length(T)
    plot3(0,0,0,'d', 'MarkerEdgeColor', 'black','MarkerFaceColor','r','MarkerSize', 10)
    plot3(err_r(i,1)/1000,err_r(i,2)/1000,err_r(i,3)/1000,'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','b','MarkerSize', 10)
    
    
    if err_r_mag(i)/1000/ERR_r >= 1
        colour = [1 0 0];
    else
        colour = [err_r_mag(i)/1000/ERR_r (1-err_r_mag(i)/1000/ERR_r)*0.5 0];
    end
    words = num2str((err_r_mag(i)/1000),'%.3f');
    text(.8*ERR_r,-ERR_r,.8*ERR_r,['Error = ' words 'km'],'Color',colour,'FontSize',18);
    text(.8*ERR_r,-ERR_r,.6*ERR_r,['Orbits = ' num2str(T(i)/perO,'%.1f')],'FontSize',18);
    
    line([0 err_r(i,1)/1000], [0 err_r(i,2)/1000], [0 err_r(i,3)/1000], 'LineStyle','-','Color', colour,'LineWidth',2 )
    
    light('Position',[1 0 0],'Style','infinite')
    pause(0.002)
    
    f = getframe(gcf);
    writeVideo(v, f);
    cla
end
close(v)


%% Orbit Animations

%Initializing the Drawing Space
set(gcf,'Renderer','zbuffer','Menubar','default','Name','Orbit Visualization', ... 
    'NumberTitle','off','Position',[300,50,650,650], ... 
    'Color',[01 01 01]); 
lim=(1+ecc)*sma;%Setting the limits of the graph
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
title('Orbital Visualization');


%Plotting the Earth

v = VideoWriter([pwd '\animations\Orbit_.mpeg']); 

open(v)

GMST_deg = radtodeg(gst);
for i=1:5:length(T)
equat_rad=rE;
    polar_rad=6356752.3142;
    [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
    load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';     
    pro.Cdata = topo2;
    earth= surface(xx,yy,zz,pro);
    colormap(topomap1) 
rotate (earth, [0 0 1], GMST_deg(i));
Xaxis= line([0 lim],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis= line([0 0],[0 lim],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
rotate (Xaxis, [0 0 1], GMST_deg(i));
rotate (Yaxis, [0 0 1], GMST_deg(i));
light('Position',r_ES(i,1:3),'Style','infinite')

line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')

plot3(r_all(1:i,1),r_all(1:i,2),r_all(1:i,3),'r')
plot3(r_all(i,1),r_all(i,2),r_all(i,3),'d', 'MarkerEdgeColor', 'black','MarkerFaceColor','r','MarkerSize', 10)
plot3(r_J2(1:i,1),r_J2(1:i,2),r_J2(1:i,3),'b')
plot3(r_J2(i,1),r_J2(i,2),r_J2(i,3),'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','b','MarkerSize', 10)
hold on


f = getframe(gcf);
writeVideo(v, f);
cla

end

close(v)


close all


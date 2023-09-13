% Author: AstroGentoo (at khu.ac.kr)
% Date: 2023.09.14.
% Purpose: Predict the equilibrium shapes of crystals based on the Wulff
% construction method in 2-D

% Description:
% The surface free energy of a crystal depends on its crystallographic
% orientation. The Wulff construction is a geometrical procedure used to
% determine the equilibrium shape of a crystal that minimizes the total
% surface energy. To obtain this shape, draw planes perpendicular to the
% normal unit vectors on the surface free energy. The inner envelope formed
% by these planes corresponds to the equilibrium crystal shape.

% Adapted from:
% Wang, L., Shirodkar, S.N., Zhang, Z. et al.
% Defining shapes of two-dimensional crystals with undefinable edge energies
% Nat Comput Sci 2, 729–735 (2022)
% https://doi.org/10.1038/s43588-022-00347-5
% Code availability
% https://zenodo.org/record/7130224
% MATLAB_codes.zip

clear all;

g=1.0; % base
amp=0.25; % amplitude of the spherical Gaussian
shp=5.0; % sharpness of the spherical Gaussian

%delta_s=0.05; % amplitude of the trigonometric function

coord=zeros(360,2);

% surface energy
for i=1:360
    theta=i*2.0*pi/360.0; % polar angle

    % cartesian coordinate
    v1=cos(theta);
    v2=sin(theta);

    % spherical Gaussian function
    g_g=g*(1.0-amp*exp(shp*(1.0*v1+0.0*v2-1.0)));
    g_g=g_g+g*(-amp*exp(shp*(0.0*v1+1.0*v2-1.0)));
    g_g=g_g+g*(-amp*exp(shp*(-1.0*v1+0.0*v2-1.0)));
    g_g=g_g+g*(-amp*exp(shp*(0.0*v1-1.0*v2-1.0)));
    g_g=g_g+g*(amp*exp(shp*(0.707*v1+0.707*v2-1.0)));
    g_g=g_g+g*(amp*exp(shp*(-0.707*v1-0.707*v2-1.0)));
    g_g=g_g+g*(amp*exp(shp*(0.707*v1-0.707*v2-1.0)));
    g_g=g_g+g*(amp*exp(shp*(-0.707*v1+0.707*v2-1.0)));

%    g_g=g*(1.0+delta_s*cos(4.0*theta)); % trigonometric function
    
    % polar coordinates
    coord(i,1)=g_g*v1;
    coord(i,2)=g_g*v2;
end

% plot of the surface energy
figure
plot(coord(:,1),coord(:,2))
axis square

% Wulff shape
w=zeros(3001,3600);
w(:,1)=-1.5:0.001:1.5; % range and interval

for i=2:360
    for j=1:3001
        % equation of normal
        w(j,i)=-(coord(i-1,1)/coord(i-1,2))*(w(j,1)-coord(i-1,1))+coord(i-1,2);
    end
end

% plot of the equilibrium shape given by Wulff construction
figure
p1=plot(coord(:,1),coord(:,2),"red",'LineWidth',2); hold on;
for k=2:360
    p2=plot(w(:,1),w(:,k),"blue");
    p2.Color(4)=0.1;
    hold on;
end
axis square
xlim([-1.5 1.5])
ylim([-1.5 1.5])
view(-90,90);

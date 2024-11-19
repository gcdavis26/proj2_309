%Geofrey Davis, Carlos Andrade: Group 107
%We have completed this assignment with integrity
clc
clear
close

%Defining ECI basis vectors
I = [1;0;0];
J = [0;1;1];
K = [0;0;1]; 

%Initial conditions
r1E = [-4357.1; 7885.8; -2225.3];
v1E = [-4.224;-5.718; -4.864];
t1 = 1564;
t2 = 32548;
mu = 3.986 * 10^5;

h = cross(r1E,v1E); %Angular momentum
n = cross(K,h)/norm(cross(K,h));

%Cap omega
capOmega = atand(dot(n,J)/dot(n,I));
if dot(n,I) < 0 %Quadrant correction
    capOmega = capOmega + 180;
end

%eccentricity
e = cross(v1E, h)/mu - r1E / norm(r1E);
e_hat = e/norm(e);

%Low omega
lowOmega = acosd(dot(n,e_hat));
if dot(e,K) < 0 %Quadrant correction
    lowOmega = -lowOmega;
end

%inclination
i = acosd(dot(K,h)/norm(h));

%semi-major axis
%E = v^2 / 2 - mu / r = -mu / 2a
%a = (-mu / 2) / (v^2 / 2 - mu / r)
a = (-mu / 2) / (norm(v1E)^2 / 2 - mu / norm(r1E));

C_EP = dcm3(lowOmega) * dcm1(i) * dcm3(capOmega);

r1P = C_EP * r1E;
v1P = C_EP * v1E;

%time of periapsis passage






function [dcm] = dcm1(angle)
dcm = [1,0,0;0,cosd(angle),sind(angle);0, -sind(angle), cosd(angle)];
end

function [dcm] = dcm2(angle)
dcm = [cosd(angle),0,-sind(angle);0,1,0; sind(angle),0,cosd(angle)];
end

function [dcm] = dcm3(angle)
dcm = [cosd(angle),sind(angle),0;-sind(angle),cosd(angle),0;0,0,1];
end


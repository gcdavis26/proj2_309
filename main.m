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
e_mag = norm(e);

%Low omega
lowOmega = acosd(dot(n,e)/e_mag);
if dot(e,K) < 0 %Quadrant correction
    lowOmega = -lowOmega;
end

%inclination
i = acosd(dot(K,h)/norm(h));

%semi-major axis
%E = v^2 / 2 - mu / r = -mu / 2a
%a = (-mu / 2) / (v^2 / 2 - mu / r)
a = (-mu / 2) / (norm(v1E)^2 / 2 - mu / norm(r1E));
p = a*(1-e_mag^2);

C_EP = dcm3(lowOmega) * dcm1(i) * dcm3(capOmega);

r1P = C_EP * r1E;
v1P = C_EP * v1E;

%time of periapsis passage
theta = acosd(1/e_mag * (p/norm(r1P) - 1));
if dot(r1P, v1P)<0
    theta = -theta;
end

E1 = 2*atan(sqrt((1-e_mag)/(1+e_mag))*tand(theta/2));
tp = t1 - sqrt(a^3/norm(r1P))* (E1 - sin(E1));
T = 2*pi*sqrt(a^3/mu);
while tp - t1 > 0 %Making sure tp is before t1 
    tp = tp - T;
end


%Newton Raphson for t2 
delta = 1 * 10^-4;
M = sqrt(mu / a^3) * (t2-tp);
f = @(x) x - e_mag*sin(x) - M;
fprime = @(x) 1-e_mag*cos(x);
Eold = M;
error = delta * 2;
while error > delta
    Enew = Eold - f(Eold)/fprime(Eold);
    error = abs(Eold - Enew);
    Eold = Enew;
end

%Finding theta2, v2, r2
E2 = Eold; 
theta2 = atand(tan(E2/2) / sqrt((1-e_mag)/(1+e_mag)))*2;

r2 = norm(h)^2 / (mu * (1+e_mag * cosd(theta2)));
r2P = [r2*cosd(theta2);r2*sind(theta2);0];
v2P = [-sqrt(mu/p) *sind(theta2); sqrt(mu/p)*(e_mag + cosd(theta2));0];

% Transforming back to perifocal
r2E = C_EP' * r2P;
v2E  = C_EP' * v2P;

fprintf("a: %.6f km\n", a)
fprintf("e: %.6f\n", e_mag)
fprintf("i: %.6f degrees\n", i)
fprintf("Ω: %.6f degrees\n", capOmega)
fprintf("ω: %.6f degrees\n", lowOmega)
fprintf("tp: %.6f seconds\n", tp)

fprintf('C_EP = [%8.6f %8.6f %8.6f ] \n', C_EP(1,:), C_EP(2,:), C_EP(3,:) )
fprintf('r2 in equitorial:  [%.6f, %.6f, %.6f] km\n', r2E)
fprintf('v2 in equitorial:  [%.6f, %.6f, %.6f] km\n', v2E)


function [dcm] = dcm1(angle)
dcm = [1,0,0;0,cosd(angle),sind(angle);0, -sind(angle), cosd(angle)];
end

function [dcm] = dcm2(angle)
dcm = [cosd(angle),0,-sind(angle);0,1,0; sind(angle),0,cosd(angle)];
end

function [dcm] = dcm3(angle)
dcm = [cosd(angle),sind(angle),0;-sind(angle),cosd(angle),0;0,0,1];
end


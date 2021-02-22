%% Housekeeping

clc;
clear;
close all;

%% Gibbs Call and Initialiazation

r1 = [-105400;157900;-152];
r2 = [-146100;108100;-226.5];
r3 = [-165200;42540;-267.3];
mu = 2.654e5; % mu of the sun

% r vectors given in AU converted to km
% r1 = [0.5887;-0.2206;0.0239]*149597900;
% r2 = [0.5027;0.2289;0.0436]*149597900;
% r3 = [0.3243;0.456;0.0453]*149597900;
% mu = 132712e6; % mu of the sun

[v2, ierr] = gibbs(r1, r2, r3, mu);

%% Find keplarian components
h = cross(r2,v2);
n = cross([0;0;1],h);
e = cross(v2,h)/mu - r2/norm(r2);
epsilon = (norm(v2)^2)/2 - mu/norm(r2);
a = - mu/(2*epsilon);
T = 2*pi*sqrt(a^3/mu);
e_norm = norm(e);
i = acosd(h(3)/norm(h));
if n(2) < 0
    Omega = 360 - acosd(n(1)/norm(n));
else
    Omega = acosd(n(1)/norm(n));
end
if e(3) < 0
    omega = 360 - acosd(dot(n,e)/(norm(n)*norm(e)));
else
    omega = acosd(dot(n,e)/(norm(n)*norm(e)));
end
if dot(r2,v2) < 0
    theta =  360 - acosd(dot(e,r2)/(norm(e)*norm(r2)));
else
    theta =  acosd(dot(e,r2)/(norm(e)*norm(r2)));
end



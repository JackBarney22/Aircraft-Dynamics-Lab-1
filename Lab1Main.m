% Jack Barney
% ASEN 3128
% 8/23/22
% Lab1Main.m
clear all, clc

%% Problem 1
%% Initial Conditions 
tspan = [0 1];
x0 = 0.1;
y0 = 0.1;
z0 = 0.1;
X0 = [x0;y0;z0];
%% Calling ODE
[t,X] = ode45(@(t,X) ODEFUN(t,X),tspan,X0);
%% Plotting ODE
% 3D Plot
% figure(1)
% hold on
% grid on
% plot3(X(:,1),X(:,2),X(:,3))
% title('3D Plot')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% hold off
% X Plot
figure(2)
subplot(3,1,1)
plot(t,X(:,1))
hold on
grid on
title('X vs Time')
xlabel('Time')
ylabel('X Value')
hold off
% Y Plot
subplot(3,1,2)
plot(t,X(:,2))
hold on
grid on
title('Y vs Time')
xlabel('Time')
ylabel('Y Value')
hold off
% Z Plot
subplot(3,1,3)
plot(t,X(:,3))
hold on
grid on
title('Z vs Time')
xlabel('Time')
ylabel('Z Value')
hold off
%% Problem 2
%% Initial Conditions
tspan = [0 10];
g = 9.81; %Meters per Senond Squared
m0 = 0.03; %Kilograms
Cd = 0.6; 
d = 0.03; %Meters
A = pi * (0.5*d)^2; % Square Meters
rho = 1.225; %kg/m^3
x0 = 0; % Meters North
y0 = 0; % Meters East
z0 = 0; % Meters Down
vx0 = 0; % m/s North
vy0 = 20; % m/s East
vz0 = -20; % m/s Down
wind = [0;0;0]; 
X0 = [x0;y0;z0;vx0;vy0;vz0];
%% Calling ODE Function
[t,X] = ode45(@(t,X) objectEOM(t,X,rho,Cd,A,m0,g,wind),tspan,X0);
%% Plotting ODE
figure(3)
hold on 
grid on
plot3(X(:,1),X(:,2),X(:,3))
plot3(X(end,1),X(end,2),X(end,3),'ro')
title('Trajectory of Golf Ball')
xlabel('North Distance(m)')
ylabel('East Distance(m)')
zlabel('Height (m)')
hold off
%% Affect of Wind
n = 30;
Xs = zeros(1,n);
Ys = zeros(1,n);
Dist = zeros(1,n);
wind_v = linspace(1,n,n);
for i = 1:n
    wind_1 = [i;0;0];
    Xs(1,i) = X(end,1);
    Ys(1,i) = X(end,2);
    Dist(1,i) = sqrt(Xs(i).^2+Ys(i).^2);
    [t,X] = ode45(@(t,X) objectEOM(t,X,rho,Cd,A,m0,g,wind_1),tspan,X0);
end
figure(4)
plot(wind_v,Dist)
hold on 
grid on
title('Final Distance vs. Wind Speed(In The North Direction)')
xlabel('Wind Speed (m/s)')
ylabel('Distance Traveled (m)')
hold off
%% Constrained Kinetic Energy
KE = 0.5*m0*(vx0^2+vy0^2+vz0^2); %Total Kinetic Energy
m_vec = linspace(0.001,0.2,200); %Vector containing different mass values
speed_vec = sqrt(2*KE./m_vec); %Vector containing speed corresponding to each mass
k = length(m_vec);
Dist_KE = zeros(1,k);
for i = 1:k
m = m_vec(i);
vy1(i) = speed_vec(i)*cos(pi/4); % m/s East
vz1(i) = -speed_vec(i)*cos(pi/4); % m/s Down 
X0 = [x0;y0;z0;vx0;vy1(i);vz1(i)];
Xs_KE(1,i) = X(end,1);
Ys_KE(1,i) = X(end,2);
Dist_KE(1,i) = sqrt(Xs_KE(i).^2+Ys_KE(i).^2);
[t,X] = ode45(@(t,X) objectEOM(t,X,rho,Cd,A,m,g,wind),tspan,X0);
end
figure(5)
plot(m_vec(1:199),Dist_KE(2:end))
hold on 
grid on
title('Distance vs. Varying Mass (Constant KE)')
xlabel('Mass (kg)')
ylabel('Distance (m)')
hold off
%% ODE Function Problem 2
function xdot = objectEOM(t,X,rho,Cd,A,m,g,wind) 
% Extract state vector
pos = [X(1);X(2);X(3)]; % Index corresponding position
v = [X(4);X(5);X(6)]; % Index corresponding to velocity
v_rel = v-wind; % Relative Wind
v_rel_mag = sqrt(v_rel(1)^2+v_rel(2)^2+v_rel(3)^2); %Magnitude of Relative Wind
h = v_rel./v_rel_mag; %Heading Vector

fGrav = [0;0;m*g];
fDrag = h.*(-1*0.5*rho*v_rel_mag^2*Cd*A);
v = [X(4);X(5);X(6)];

if pos(3) > 0 %Ending Condition
     fgrav = [0;0;0];
     fdrag = [0;0;0];
     v = [0;0;0];
end 

fnet = fGrav + fDrag;
a = fnet./m;
xdot = [v(1);v(2);v(3);a(1);a(2);a(3)];
end
%% ODE Function Problem 1
function dX = ODEFUN(t,X)

xd = X(1) + 2*X(2) + X(3);
yd = X(1) - 5*X(3);
zd = X(1)*X(2) - X(2).^2 +3*X(3).^3;
dX = [xd;yd;zd];

end

% Jack Barney
% ASEN 3128
% 8/23/22
% Lab1Problem1.m
clear all, clc

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
figure(1)
hold on
grid on
plot3(X(:,1),X(:,2),X(:,3))
title('3D Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off
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
%% ODE Function
function dX = ODEFUN(t,X)

xd = X(1) + 2*X(2) + X(3);
yd = X(1) - 5*X(3);
zd = X(1)*X(2) - X(2).^2 +3*X(3).^3;
dX = [xd;yd;zd];

end

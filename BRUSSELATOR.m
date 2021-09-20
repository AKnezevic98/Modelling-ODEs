%% THE BRUSSELATOR
%{
This following system of equations model the changes in concentration of chemicals X and Y which are in a chemical reaction:
dX/dt = 1 - (b+1)*X + a*X^2*Y
dY/dt = b*X - a*^2*Y
This code is a numerical integration of these equations.
%}

%% THE PARAMETERS
X0 = 0.9;
Y0 = 2.1;
a = 1;
b = 2;
T = 20;
tspan = [0, T];
xs = linspace(0,T,1000);

%% SOLVER OPTIONS
options = odeset('RelTol', 1e-6, 'AbsTol', 10e-6, 'Vectorized', 'on');
FS = 15; % FONT SIZE
LW = 2; % LINE WIDTH

%% INTEGRATION

f = figure();

% INTEGRATION OF THE SYSTEM
subplot(2,2,1);
y0 = [X0, Y0];
sol = ode45(@(t,y) odefunc(t,y,a,b), tspan, y0, options); ys = deval(sol, xs);
plot(xs, ys(1,:), 'LineWidth', LW);
hold on;
plot(xs, ys(2,:), 'LineWidth', LW);
tt = (1 - (b-1)/b);
T0 = 2*asin(2*tt)/(a*tt);
ymax = ylim;
plot([T0, T0],[ymax(1), ymax(2)],  'LineWidth', LW, 'LineStyle', '--');
plot([2*T0, 2*T0],[ymax(1), ymax(2)],  'LineWidth', LW, 'LineStyle', '--');
title('a)', 'FontSize', FS);
legend({'Chemical X', 'Chemical Y', sprintf('T = %.2f', T)});
grid on;
xlabel('Time t');
ylabel('Concentration');
ax = gca;
ax.FontSize = FS;

% NULLCLINES
subplot(2,2,2);
xmin = 0.8; xmax = 1.2; ymin = b/a-0.2; ymax = b/a+0.2;
xs = linspace(xmin, xmax, 100);
ys = ((b+1)*xs-1)./(a*xs.^2);
plot(xs, ys, 'LineWidth', LW);
hold on;
ys = b./(a*xs);
plot(xs, ys, 'LineWidth', LW);
scatter([1],[b/a],'Filled');
y0 = [0.9, b/a+0.1];
scatter([y0(1)], [y0(2)], 'Filled');
ts = linspace(0,100*T,100000); tspan = [0, 100*T];
sol = ode45(@(t,y) odefunc(t,y,a,b), tspan, y0, options); ys = deval(sol, ts);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);

xs = linspace(xmin, xmax, 10);
ys0 = linspace(ymin, ymax, 10);
[Xs,Ys] = meshgrid(xs,ys0);
Us = 1 - (b+1)*Xs + a*Xs.^2.*Ys;
Vs = b*Xs - a*Xs.^2.*Ys;
quiver(xs, ys0, Us, Vs, 'AutoScaleFactor', 1);

ylim([ymin, ymax]);
xlim([xmin, xmax]);
title('b)', 'FontSize', FS);
legend({'X nullcline', 'Y nullcline', sprintf('Fixed point (1, %.2f)', b/a), 'Starting point'});
grid on;
xlabel('Concentration of X');
ylabel('Concentration of Y');
ax = gca;
ax.FontSize = FS;


% BIFURCATION
subplot(2,2,3);
b = 1 + a + 0.1;
xmin = 0.5; xmax = 1.5; ymin = b/a-0.5; ymax = b/a+0.5;
xs = linspace(xmin, xmax, 100);
ys = ((b+1)*xs-1)./(a*xs.^2);
plot(xs, ys, 'LineWidth', LW);
hold on;
ys = b./(a*xs);
plot(xs, ys, 'LineWidth', LW);
scatter([1],[b/a],'Filled');
y0 = [0.9, b/a+0.1]; 
scatter([y0(1)], [y0(2)], 'Filled');
ts = linspace(0,100*T,100000); tspan = [0, 100*T];
sol = ode45(@(t,y) odefunc(t,y,a,b), tspan, y0, options); ys = deval(sol, ts);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);

xs = linspace(xmin, xmax, 10);
ys0 = linspace(ymin, ymax, 10);
[Xs,Ys] = meshgrid(xs,ys0);
Us = 1 - (b+1)*Xs + a*Xs.^2.*Ys;
Vs = b*Xs - a*Xs.^2.*Ys;
quiver(xs, ys0, Us, Vs, 'AutoScaleFactor', 1);

ylim([ymin, ymax]);
xlim([xmin, xmax]);
title('c)', 'FontSize', FS);
legend({'X nullcline', 'Y nullcline', sprintf('Fixed point (1, %.2f)', b/a), 'Starting point'});
grid on;
xlabel('Concentration of X');
ylabel('Concentration of Y');
ax = gca;
ax.FontSize = FS;

subplot(2,2,4);
b = 1 + a - 0.1;
xmin = 0.5; xmax = 1.5; ymin = b/a-0.5; ymax = b/a+0.5;
xs = linspace(xmin, xmax, 100);
ys = ((b+1)*xs-1)./(a*xs.^2);
plot(xs, ys, 'LineWidth', LW);
hold on;
ys = b./(a*xs);
plot(xs, ys, 'LineWidth', LW);
scatter([1],[b/a],'Filled');
y0 = [0.9, b/a+0.1]; 
scatter([y0(1)], [y0(2)], 'Filled');
ts = linspace(0,100*T,10000); tspan = [0, 100*T];
sol = ode45(@(t,y) odefunc(t,y,a,b), tspan, y0, options); ys = deval(sol, ts);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);

xs = linspace(xmin, xmax, 10);
ys0 = linspace(ymin, ymax, 10);
[Xs,Ys] = meshgrid(xs,ys0);
Us = 1 - (b+1)*Xs + a*Xs.^2.*Ys;
Vs = b*Xs - a*Xs.^2.*Ys;
quiver(xs, ys0, Us, Vs, 'AutoScaleFactor', 1);

ylim([ymin, ymax]);
xlim([xmin, xmax]);
title('d)', 'FontSize', FS);
legend({'X nullcline', 'Y nullcline', sprintf('Fixed point (1, %.2f)', b/a), 'Starting point'});
grid on;
xlabel('Concentration of X');
ylabel('Concentration of Y');
ax = gca;
ax.FontSize = FS;


%% THE SYSTEMS OF EQUATIONS
function dy = odefunc(t,y,a,b)
dy = [0; 0];
dy(1) = 1 - (b+1)*y(1) + a*y(1).^2.*y(2);
dy(2) = b*y(1) - a*y(1).^2.*y(2);
end

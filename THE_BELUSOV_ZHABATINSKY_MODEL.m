%% THE BELOUSOV ZHABOTINSKY MODEL
%{
This following system of equations model the changes in concentration of chemicals X, Y and Z which are in a chemical reaction:
dX/dt = a * (Y + X*(1 - b*X - Y))
dY/dt = (Z - (1 + X)*Y)/a
dZ/dt = g*(X - Z)
This code is a numerical integration of these equations.
%}

%% THE PARAMETERS
X0 = 1;
Y0 = 2;
Z0 = 3;
alpha = 77.27;
beta = 8.375e-6;
gamma = 1.161;
T = 50;
tspan = [0, T];

%% SOLVER OPTIONS
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Vectorized', 'on', 'Stats', 'on', 'Jacobian', @(t,y) odejac(t,y,alpha,beta,gamma));
FS = 15; % FONT SIZE
LW = 2; % LINE WIDTH

%% INTEGRATION

f = figure();

% INTEGRATION USING ODE45
subplot(2,3,1);
y0 = [X0, Y0, Z0];
disp('***ODE45:');
sol = ode45(@(t,y) odefunc(t,y,alpha,beta,gamma), tspan, y0, options);
%ys = deval(sol, xs); xs = linspace(0,T,1000);
ys = sol.y; xs = sol.x;
dxs = xs(2:end) - xs(1:end-1);
dxmax = max(dxs);
p1 = plot(xs, dxmax*ys(1,:)/max(ys(1,:)), '-', 'LineWidth', LW);
hold on;
p2 = plot(xs, dxmax*ys(2,:)/max(ys(2,:)), '-', 'LineWidth', LW);
p3 = plot(xs, dxmax*ys(3,:)/max(ys(3,:)), '-', 'LineWidth', LW);
p4 = plot(xs(1:end-1), dxs, 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4)=0.5; p4.Color(4) = 0.5; % OPACITY
title('a)', 'FontSize', FS);
legend({'Chemical X', 'Chemical Y', 'Chemical Z', 'Step Size'});
grid on;
xlabel('Time t');
ylabel('Concentration');
ax = gca;
ax.FontSize = FS;

% INTEGRATION USING ODE113
subplot(2,3,2);
disp('***ODE113:');
sol = ode113(@(t,y) odefunc(t,y,alpha,beta,gamma), tspan, y0, options);
%ys = deval(sol, xs); xs = linspace(0,T,1000);
ys = sol.y; xs = sol.x;
dxs = xs(2:end) - xs(1:end-1);
dxmax = max(dxs);
p1 = plot(xs, dxmax*ys(1,:)/max(ys(1,:)), '-', 'LineWidth', LW);
hold on;
p2 = plot(xs, dxmax*ys(2,:)/max(ys(2,:)), '-', 'LineWidth', LW);
p3 = plot(xs, dxmax*ys(3,:)/max(ys(3,:)), '-', 'LineWidth', LW);
p4 = plot(xs(1:end-1), dxs, 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4)=0.5; p4.Color(4) = 0.5; % OPACITY
title('b)', 'FontSize', FS);
legend({'Chemical X', 'Chemical Y', 'Chemical Z', 'Step Size'});
grid on;
xlabel('Time t');
ylabel('Concentration');
ax = gca;
ax.FontSize = FS;

% INTEGRATION USING ODE15s
subplot(2,3,3);
disp('***ODE15s:');
sol = ode15s(@(t,y) odefunc(t,y,alpha,beta,gamma), tspan, y0, options);
%ys = deval(sol, xs); xs = linspace(0,T,1000);
ys = sol.y; xs = sol.x;
dxs = xs(2:end) - xs(1:end-1);
dxmax = max(dxs);
p1 = plot(xs, dxmax*ys(1,:)/max(ys(1,:)), '-', 'LineWidth', LW);
hold on;
p2 = plot(xs, dxmax*ys(2,:)/max(ys(2,:)), '-', 'LineWidth', LW);
p3 = plot(xs, dxmax*ys(3,:)/max(ys(3,:)), '-', 'LineWidth', LW);
p4 = plot(xs(1:end-1), dxs, 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4)=0.5; p4.Color(4) = 0.5; % OPACITY
title('c)', 'FontSize', FS);
legend({'Chemical X', 'Chemical Y', 'Chemical Z', 'Step Size'});
grid on;
xlabel('Time t');
ylabel('Concentration');
ax = gca;
ax.FontSize = FS;

% 3D PLOT
subplot(2,3,5);
plot3(ys(1,:), ys(2,:), ys(3,:), 'LineWidth', LW);
hold on;
scatter3([X0], [Y0], [Z0],'Filled');
title('e)', 'FontSize', FS);
grid on;
xlabel('Concentration of X');
ylabel('Concentration of Y');
zlabel('Concentration of Z');
ax = gca;
ax.FontSize = FS;

% JACOBIAN COND
subplot(2,3,6);
xlen = length(xs);
Conds = zeros(xlen);
for i = 1:xlen
    Jac = odejac(xs(i), ys(:,i), alpha, beta, gamma);
    Conds(i) = cond(Jac);
end
plot(xs, Conds, 'LineWidth', LW);
title('f)', 'FontSize', FS);
grid on;
xlabel('Time t');
ylabel('Jacobian condition');
ax = gca;
ax.FontSize = FS;

% INTEGRATION USING ODE23s
subplot(2,3,4);
disp('***ODE23s:');
sol = ode23s(@(t,y) odefunc(t,y,alpha,beta,gamma), tspan, y0, options);
%ys = deval(sol, xs); xs = linspace(0,T,1000);
ys = sol.y; xs = sol.x;
dxs = xs(2:end) - xs(1:end-1);
dxmax = max(dxs);
p1 = plot(xs, dxmax*ys(1,:)/max(ys(1,:)), '-', 'LineWidth', LW);
hold on;
p2 = plot(xs, dxmax*ys(2,:)/max(ys(2,:)), '-', 'LineWidth', LW);
p3 = plot(xs, dxmax*ys(3,:)/max(ys(3,:)), '-', 'LineWidth', LW);
p4 = plot(xs(1:end-1), dxs, 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4)=0.5; p4.Color(4) = 0.5; % OPACITY
title('d)', 'FontSize', FS);
legend({'Chemical X', 'Chemical Y', 'Chemical Z', 'Step Size'});
grid on;
xlabel('Time t');
ylabel('Concentration');
ax = gca;
ax.FontSize = FS;


%% THE SYSTEMS OF EQUATIONS AND THE JACOBIAN
function dy = odefunc(t,y,a,b,g)
dy = [0; 0; 0];
dy(1) = a*(y(2) + y(1).*(1 - b*y(1) - y(2)));
dy(2) = (y(3) - (1+y(1)).*y(2))/a;
dy(3) = g*(y(1)-y(3));
end

function J = odejac(t,y,a,b,g)
J = zeros([3,3]);
J(1,1) = a*(1-2*b*y(1) - y(2));
J(1,2) = a*(1-y(1));
J(2,1) = -y(2)/a;
J(2,2) = -(1+y(1))/a;
J(2,3) = 1/a;
J(3,1) = g;
J(3,3) = -g;
end
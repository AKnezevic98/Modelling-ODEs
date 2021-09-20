%% THE BELOUSOV ZHABOTINSKY MODEL
%{
This following system of equations is a simplified weather model introduced by Lorenz in 1963:
dX/dt = a * (Y + X*(1 - b*X - Y))
dY/dt = (Z - (1 + X)*Y)/a
dZ/dt = g*(X - Z)
X is the rate of convection, Y is horizontal temperature variation and Z is vertical temperature variation.
This code is a numerical integration of these equations.
%}

%% THE PARAMETERS
X0 = -8;
Y0 = -8;
Z0 = 27;
beta = 8/3;
rho = 50;
sigma = 10;

%% CALCULATING THE TIME TO INTEGRATE
t0 = [-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1];
t1 = [sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), rho-1]; % FIXED POINTS
y0 = [X0, Y0, Z0];
if norm(y0-t0)<=norm(y0-t1)
    J = odejac(0,t0, beta, rho, sigma);
else
    J = odejac(0,t1, beta, rho, sigma);
end
freqs = imag(eig(J)); % FREQUENCIES
freq = max(freqs); % LARGEST FREQUENCY

T = 2*pi/freq; % CHARACTERISTIC TIME
N = 10; % NUMBER OF LOOPS I WANT
tend = 2*10;
tspan = [0, tend];

%% SOLVER OPTIONS
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Vectorized', 'on', 'Jacobian', @(t,y) odejac(t,y,beta,rho,sigma));
FS = 15; % FONT SIZE
LW = 2; % LINE WIDTH

%% INTEGRATION

f = figure();

% PLOTTING DIFFERENT ODEs TOL=1e-6
subplot(2,3,1);
sol = ode45(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p1 = plot(xs, ys(1,:), 'LineWidth', LW);
hold on;
sol = ode113(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p2 = plot(xs, ys(1,:), 'LineWidth', LW);
sol = ode15s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p3 = plot(xs, ys(1,:), 'LineWidth', LW);
sol = ode23s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p4 = plot(xs, ys(1,:), 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4) = 0.5; p4.Color(4) = 0.5;
legend({'ode45 X(t)', 'ode113 X(t)', 'ode15s X(t)', 'ode23s X(t)'});
title('a)', 'FontSize', FS);
grid on;
xlabel('Time t');
ylabel('Convection X');
ax = gca;
ax.FontSize = FS;

% PLOTTING DIFFERENT ODEs TOL=1e-14
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14, 'Vectorized', 'on', 'Jacobian', @(t,y) odejac(t,y,beta,rho,sigma));
subplot(2,3,2);
sol = ode45(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p1 = plot(xs, ys(1,:), 'LineWidth', LW);
hold on;
sol = ode113(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p2 = plot(xs, ys(1,:), 'LineWidth', LW);
sol = ode15s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p3 = plot(xs, ys(1,:), 'LineWidth', LW);
sol = ode23s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
p4 = plot(xs, ys(1,:), 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; p3.Color(4) = 0.5; p4.Color(4) = 0.5;
legend({'ode45 X(t)', 'ode113 X(t)', 'ode15s X(t)', 'ode23s X(t)'});
title('b)', 'FontSize', FS);
grid on;
xlabel('Time t');
ylabel('Convection X');
ax = gca;
ax.FontSize = FS;

% PLOTTING SOLUTION AND ANOTHER ONE JUST A LITTLE OFF
subplot(2,3,3);
p1 = plot3(ys(1,:), ys(2,:), ys(3,:), 'LineWidth', LW);
hold on;
scatter3([ys(1,end)], [ys(2,end)], [ys(3,end)], 'Filled');
sol2 = ode23s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, [X0+1e-15,Y0,Z0], options); xs2 = sol2.x; ys2 = sol2.y;
p2 = plot3(ys2(1,:), ys2(2,:), ys2(3,:), 'LineWidth', LW);
scatter3([ys2(1,end)], [ys2(2,end)], [ys2(3,end)], 'Filled');
scatter3([X0], [Y0], [Z0], 'Filled');
p1.Color(4) = 0.5; p2.Color(4) = 0.5;
legend({'Solution 1', '1st End Point', 'Solution 2', '2nd End Point'});
title('c)', 'FontSize', FS);
grid on;
xlabel('Convection X');
ylabel('Temperature variation Y');
zlabel('Temperature variation Z');
ax = gca;
ax.FontSize = FS;


% PLOTTING DIFFERENCES OF SOLUTIONS
subplot(2,3,4);
sol2 = ode23s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, [X0+1e-15,Y0,Z0], options);
sol = ode23s(@(t,y) odefunc(t,y,beta,rho,sigma), tspan, y0, options);
xs = sol.x; ys = sol.y;
Nmax = 10000;
xs0 = linspace(0, tend, Nmax);
ys1 = deval(sol, xs0); ys2 = deval(sol2, xs0);
dys1 = ys1(1,:) - ys2(1,:); dys2 = ys1(2,:) - ys2(2,:); dys3 = ys1(3,:) - ys2(3,:);
dy = zeros([1,Nmax]);
for i=1:Nmax
    dy(i) = norm([dys1(i), dys2(i), dys3(i)]);
end
linfit = fit(xs0', dy', 'exp1');
disp(linfit);
p1 = plot(xs0, dy, 'LineWidth', LW);
hold on;
p2 = plot(xs0, linfit.a*exp(linfit.b*xs0), 'LineWidth', LW);
legend({'Differences', 'Fit'});
p1.Color(4) = 0.5; p2.Color(4) = 0.5;
title('d)', 'FontSize', FS);
grid on;
xlabel('Time t');
ylabel('Value');
ax = gca;
ax.FontSize = FS;

% PEAKS
subplot(2,3,5);
zs = ys(3,:);
[pks,locs] = findpeaks(zs);
spks = pks(2:end);
peaksfit = fit(pks(1:end-1)', spks', 'pchipinterp');
scatter(pks(1:end-1), spks, 'Filled');
hold on;
p1 = plot(peaksfit);
p1.LineWidth = LW;
title('e)', 'FontSize', FS);
grid on;
xlabel('N-th peak');
ylabel('N+1-st peak');
ax = gca;
ax.FontSize = FS;

% PREDICTING PEAKS
subplot(2,3,6);
Npeaks = length(pks(1:end-1)); simpks = zeros([1, Npeaks]);
simpks(1) = peaksfit(pks(1));
for i=2:Npeaks
    simpks(i) = peaksfit(simpks(i-1));
end
p1 = plot(xs(locs(2:end)), pks(2:end), 'LineWidth', LW);
hold on;
p2 = plot(xs(locs(2:end)), simpks, 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5;
legend({'Simulation', 'Prediction'});
title('f)', 'FontSize', FS);
grid on;
xlabel('Time t');
ylabel('Peaks');
ax = gca;
ax.FontSize = FS;


%% THE SYSTEMS OF EQUATIONS AND THE JACOBIAN
function dy = odefunc(t,y,b,r,s)
dy = [0; 0; 0];
dy(1) = s*(y(2)-y(1));
dy(2) = y(1).*(r-y(3)) - y(2);
dy(3) = y(1).*y(2) - b*y(3);
end

function J = odejac(t,y,b,r,s)
J = zeros([3,3]);
J(1,1) = -s;
J(1,2) = s;
J(2,1) = r-y(3);
J(2,2) = -1;
J(2,3) = -y(1);
J(3,1) = y(2);
J(3,2) = y(1);
J(3,3) = -b;
end
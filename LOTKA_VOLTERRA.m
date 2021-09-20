%% LOTKA-VOLTERRA EQUATIONS
%{
This following system of equations model the changes in populations H of pray and C of predators:
dH/dt = a*H - b*H*C
dC/dt = c*H*C - d*C
This code is a numerical integration of these equations.
%}

%% THE PARAMETERS
H0 = 600;
C0 = 400;
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
T = 20;
tspan = [0, T];
xs = linspace(0,T,1000);

%% SOLVER OPTIONS
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4, 'Vectorized', 'on');
FS = 15; % FONT SIZE
LW = 2; % LINE WIDTH

%% INTEGRATION

f = figure();

% FIRST INTEGRATION FOR H(t) AND C(t)
subplot(2,3,1);
y0 = [H0, C0]; sol = ode45(@(t,y) odefunc1(t,y,a,b,c,d), tspan, y0, options); ys = deval(sol, xs);
plot(xs, ys(1,:), 'LineWidth', LW);
hold on;
plot(xs, ys(2,:), 'LineWidth', LW);
title('a)', 'FontSize', FS);
grid on;
xlabel('Time');
ylabel('Number of individuals');
ax = gca;
ax.FontSize = FS;

% STOCHASTIC
ts = [0,];
hs = [H0,]; cs = [C0,];
while ts(end)<T
    r = [a*hs(end), b*hs(end)*cs(end), d*cs(end)]; %RATES
    rt = sum(r); %TOTAL RATE
    ts(end+1) = ts(end) - log(random('Uniform',0,1))/rt;%NEXT TIME STEP
    rn = random('Uniform',0,1);
    if rn < r(1)/rt
        hs(end+1) = hs(end)+1;
        cs(end+1) = cs(end);
    elseif rn < sum(r(1:2))/rt
        hs(end+1) = hs(end)-1;
        cs(end+1) = cs(end)+1;
    else
        hs(end+1) = hs(end);
        cs(end+1) = cs(end)-1;
    end
end
plot(ts(1:end-1),hs(1:end-1), 'LineWidth', LW);
plot(ts(1:end-1),cs(1:end-1), 'LineWidth', LW);
legend({'Deterministic H', 'Deterministic C', 'Stochastic H', 'Stochastic C'});



% SEVERAL INTEGRATIONS FOR PHASE SPACE PLOTS
subplot(2,3,2);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
hold on;
plot(hs(1:end-1),cs(1:end-1), 'LineWidth', LW);
H01 = 600; C01 = 600; y0 = [H01, C01]; sol = ode45(@(t,y) odefunc1(t,y,a,b,c,d), tspan, y0, options); ys = deval(sol, xs);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
H02 = 400; C02 = 600; y0 = [H02, C02]; sol = ode45(@(t,y) odefunc1(t,y,a,b,c,d), tspan, y0, options); ys = deval(sol, xs);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
H03 = ceil(d/c); C03 = ceil(a/b);
scatter(H03, C03, 'filled');
title('b)', 'FontSize', FS);
legend({sprintf('H_0,C_0=%d,%d',H0,C0), sprintf('Stoch. %d,%d',H0,C0), sprintf('H_0,C_0=%d,%d',H01,C01), sprintf('H_0,C_0=%d,%d',H02,C02), sprintf('H_0,C_0=%d,%d',H03,C03)});
grid on;
xlabel('Pray population H');
ylabel('Predator population C');
ax = gca;
ax.FontSize = FS;

% INTEGRATIONS USING THE DIMENSIONLESS EQUATIONS
subplot(2,3,3);
hold on;
alpha1 = 9/4; h0 = 2/3; c0 = 1; y0 = [h0, c0]; sol = ode45(@(t,y) odefunc2(t,y,alpha1), tspan, y0, options); ys = deval(sol, xs);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
alpha2 = 4/4; sol = ode45(@(t,y) odefunc2(t,y,alpha2), tspan, y0, options); ys = deval(sol, xs);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
alpha3 = 16/4; sol = ode45(@(t,y) odefunc2(t,y,alpha3), tspan, y0, options); ys = deval(sol, xs);
plot(ys(1,:), ys(2,:), 'LineWidth', LW);
scatter(1, 1, 'filled');
title('c)', 'FontSize', FS);
grid on;
legend({sprintf('\\alpha=%.2f',alpha1), sprintf('\\alpha=%.2f',alpha2), sprintf('\\alpha=%.2f',alpha3), sprintf('H_0,C_0=%i,%i',1,1)});
xlabel('Pray population H');
ylabel('Predator population C');
ax = gca;
ax.FontSize = FS;

% INTEGRATIONS WITH ADDED CONDITIONS
subplot(2,3,4);
hold on;
c = 10; s = 7; T = 1000; tspan = [0, T];
alpha = 9/4; h01 = 1; c01 = 1; y0 = [h01, c01]; sol = ode45(@(t,y) odefunc3(t,y,alpha,c,s), tspan, y0, options); xs = linspace(0,T,100000); ys = deval(sol, xs);
p1 = plot(ys(1,:), ys(2,:), 'LineWidth', LW);
h02 = 0.2; c02 = 1; y0 = [h02, c02]; sol = ode45(@(t,y) odefunc3(t,y,alpha,c,s), tspan, y0, options); ys = deval(sol, xs);
p2 = plot(ys(1,:), ys(2,:), 'LineWidth', LW);
p1.Color(4) = 0.5; p2.Color(4) = 0.5; % OPACITY
scatter([1,0.2],[1,1],'filled');
title('d)', 'FontSize', FS);
grid on;
legend({sprintf('H_0,C_0=%.2f,%.2f',h01,c01), sprintf('H_0,C_0=%.2f,%.2f',h02,c02)});
xlabel('Pray population H');
ylabel('Predator population C');
ax = gca;
ax.FontSize = FS;

% INTEGRATION OF 3-SPECIES SYSTEM
subplot(2,3,5);
hold off;
T = 100; tspan = [0, T];
alpha = 1; beta = 2; gamma = 3; delta = 5;
h0 = 1; o0 = 0.5; c0 = 0.3;
y0 = [h0, o0, c0]; sol = ode45(@(t,y) odefunc4(t,y,alpha,beta,gamma,delta), tspan, y0, options); xs = linspace(0,T,10000); ys = deval(sol, xs);
plot3(ys(1,:),ys(2,:),ys(3,:), 'LineWidth', LW);
hold on;
title('e)', 'FontSize', FS);
scatter3(h0,o0,c0,'filled');
grid on;
xlabel('Pray population H');
ylabel('Intermediary population O');
zlabel('Predator population C')
ax = gca;
ax.FontSize = FS;

% STOCHASTIC HEAT MAP
subplot(2,3,6);
hold off;
Xmax = 100; Xstep = floor(1);
X = 0:Xstep:Xmax; Y = 0:Xstep:Xmax;
Nmax = length(X); Z = zeros([Nmax,Nmax]);
Hs0 = zeros([Nmax, Nmax]); Cs0 = zeros([Nmax, Nmax]);
Hs = zeros([Nmax, Nmax]); Cs = zeros([Nmax, Nmax]); One = ones([Nmax,Nmax]);
Ntot = 100; % NUMBER OF SIMULATIONS FOR EACH STARTING POSITION
Nstep = 100; % NUMBER OF STEPS IN EACH SIMULATION
for n=1:Nmax
    for m = 1:Nmax
        Hs0(n,m) = X(n);
        Cs0(n,m) = Y(m);
    end
end
a = 0.4; b = 0.001; c = 0.001; d = 0.9;
for x = 1:Ntot
    Hs = Hs0; Cs = Cs0; % SET TO INITIAL
    for y = 1:Nstep
        r1 = a*Hs;
        r2 = b*Hs.*Cs;
        r3 = d*Cs; % RATES
        
        rt = r1 + r2 + r3; % TOTAL RATES;
        
        rn = random('Uniform',0,1,[Nmax,Nmax]); % RANDOM NUMBERS
        
        mask1 = rn < r1./rt;
        mask2 = rn < (r1+r2)./rt - mask1;
        mask3 = One - mask1 - mask2;
        mask4 = Hs == 0;
        mask5 = Cs == 0; % MASKS FOR CONDITIONS
        
        Hs = (Hs + mask1 - mask2).*(1-mask4);
        Cs = (Cs + mask2 - mask3).*(1-mask5); % EVOLUTION
    end
    mask = Hs==0;
    mask = mask + Cs==0;
    
    Z = Z + mask/Ntot;
end
[S,K] = contourf(X,Y,Z);
title('f)', 'FontSize', FS);
colorbar;
K.LineStyle = 'none';
xlabel('Initial Pray population H_0');
ylabel('Initial Predator population C_0');
ax = gca;
ax.FontSize = FS;

%cleanfigure; matlab2tikz('LOTKA-VOLTERRA_PLOT.tex'); % GET TIKZ CODE

%% THE SYSTEMS OF EQUATIONS
% NON-DIMENSIONLESS
function dy = odefunc1(t,y,a,b,c,d)
dy = [0; 0];
dy(1) = a*y(1) - b*y(1).*y(2);
dy(2) = c*y(1).*y(2) - d*y(2);
end

% DIMENSIONLESS
function dy = odefunc2(t,y,a)
dy = [0; 0];
dy(1) = y(1).*(1-y(2));
dy(2) = a*y(2).*(y(1)-1);
end

% WITH ADDED CONDITIONS
function dy = odefunc3(t,y,a,c,s)
dy = [0; 0];
dy(1) = y(1).*((1-y(1)/c)-y(2)./(1+y(1)/s));
dy(2) = a*y(2).*(y(1)./(1+y(1)/s)-1);
end

% 3 SPECIES
function dy = odefunc4(t,y,a,b,g,d)
dy = [0; 0; 0];
dy(1) = y(1).*(1-y(2)-y(3));
dy(2) = y(2).*(g*y(1) - d*y(3));
dy(3) = a*y(3).*(y(1) + b*y(2)-1);
end
% randomWalk2d.m
% Ben Zager
% Math 87 Final Project
% Simulation of 2D diffusion by random walk

%% Constants (Can change)

N = 1.0*10^(5); % number of particles
T = 300; % temperature (Kelvin)
r = 1; % radius of particles (nm)
time = 1; % time interval (s)
dt = 0.001; % length of each time step (s)

%% Constants (DO NOT CHANGE)

kB = 8.617 * 10^(-5); % Boltzmann's constant (eV/K)
mu = 5.714 * 10^(-12); % viscosity of water  eV*s/(nm^3)
D = kB*T/(4*pi*mu*r); % Stokes-Einstein Relation (nm^2/s)
L = sqrt(4*D*dt); % length of each step (nm)

%% Boundary conditions
%'U'-unbounded 'S'- square 'C'- channel
shape = 'S';
bound = 3*10^(4); % general size parameter
wid = 0.1; % width of 'L' and channel (fraction of bound)
len = 0.2; % length of channel (fraction of bound)

%% Initial position vector

lin = linspace(-bound,bound,N);

%init = zeros(N,2); % point source
init = [lin; zeros(1,N)]'; % horizontal line source through middle
%init = [bound*ones(1,N); lin]'; % vertical line source at edge

%% Potential field

a = 0; % strength of drift potential in x direction
k = [0.0 0.0]; % strength of harmonic potential
U = 1; %resting membrane potential

%% Simulation
close all;

%distr = 'bin'; % binary {-1,1}
distr = 'uni'; % uniform [-1,1]
%distr = 'norm'; % std. normal <x>=0,sd=1

if (strcmp(shape,'S')) % square
    pos = square(N,L,time,dt,bound,distr,init,k);
elseif (strcmp(shape,'C')) % channel
    [pos,V] = channel(N,L,time,dt,bound,wid,len,distr,init,a,U);
elseif (strcmp(shape,'U')) % unbounded
    [pos] = unbounded(N,L,time,dt,distr,init,a,k);
elseif (strcmp(shape,'A')) % absorbing
    pos = absorbing(N,L,time,dt,bound,distr,init,a,k);
end

%% Post-Simulation Plotting

t = linspace(0,time,time/dt);
figure;
plot(t,V);
title('Potential vs. Time for Active Diffusion', 'interpreter','latex');
xlabel('time (s)','interpreter','latex');
ylabel('Potential','interpreter','latex');
hold on;
plot(t,0.5*exp(-10*t));
daspect([1 1 1]);
legend('Potential','0.5e^{-10t}','interpreter','latex');
axis([0 0.7 0 0.5])

%% Plot for figure 6

potentials = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.2];
times = [10.66 2.42 1.29 0.93 0.63 0.53 0.45 0.4 0.35 0.31 0.28 0.24];

figure; 
plot(potentials,times);
title('Time to reach equilibrium vs. starting potential', 'interpreter','latex');
xlabel('Potential','interpreter','latex');
ylabel('time (s)','interpreter','latex');

%% Bivariate gaussian plot

% figure;
% mu = [0 0]; % mean for X,Y 
% Sigma = [.25 .3; .3 1]; % covariance matrix
% x = -3:.2:3; y = -3:.2:3; 
% [X,Y] = meshgrid(x,y);
% F = mvnpdf([X(:) Y(:)],mu,Sigma); 
% F = reshape(F,length(Y),length(X));
% surf(X,Y,F);
% caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis([-2 2 -2 2 0 .4])
% xlabel('x'); ylabel('y'); zlabel('Probability Density','Interpreter','latex');
% titleText = 'Bivariate Gaussian Distribution';
% title(titleText, 'Interpreter','latex')


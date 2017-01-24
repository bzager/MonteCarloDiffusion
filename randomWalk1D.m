% randomWalk1D.m 
% Ben Zager
% Math 87 Final Project
% Simulation of 1D diffusion by random walk

%% Constants

N = 1.0*10^(6); % number of particles 
kB = 8.617 * 10^(-5); % Boltzmann's constant (eV/K)
mu = 5.7138 * 10^(-12); % viscosity of water  eV*s/(nm^3)
T = 200; % temperature (Kelvin)
r = 100; % radius of particles (nm)
D = kB*T/(2*pi*mu*r); % Stokes-Einstein Relation
time = 300; % time interval (seconds)
L = sqrt(2*D); % length of each step

%% Boundary conditions

bound = 0.1*10^(5); % (nm)

%% Concentration Calculation

cells = 50;
concMeas = zeros(1,cells);
x = linspace(-bound,bound,cells);
dx = 2*bound/cells;

%% Initial conditions

pos = zeros(N,1) + 0.0*bound; % point source with some shift

%% Simulation
close all;
figure;

for t = 1:time
    step = (2 * randi([0,1],N,1)) - 1; % binary (-1 or 1)
    %step = (2 * rand(N,1)) - 1;        % uniform (-1 to 1)
    %step = randn(N,1);                  % std. normal (-1 to 1)
    pos = pos + L*step;
    pos(pos >= bound) = bound;    % boundary check for square boundary
    pos(pos <= -bound) = -bound;
    
    for i = 1:cells
        count = sum(pos > x(i) & pos < x(i)+dx);
        concMeas(i) = count/N;
    end
    
    plot(x,concMeas);
    axis([-bound,bound,0,0.4]);
    label = sprintf('Concentration at t = %d s',t);
    title(label);
    drawnow
end

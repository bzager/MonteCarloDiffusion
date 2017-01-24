%
%
%

%% Constants

N = 1*10^(5); % number of particles 
kB = 8.617 * 10^(-5); % Boltzmann's constant (eV/K)
mu = 5.7138 * 10^(-12); % viscosity of water  eV*s/(nm^3)
T = 300; % temperature (Kelvin)
r = 10; % radius of particles (nm)
D = kB*T/(6*pi*mu*r); % Einstein's diffusion relation
time = 400; % time interval (seconds)
L = sqrt(6*D); % length of each step

%%

bound = 1*10^(6); % distance of square boundaries from origin (nm)


%% Initial position vector

pos = zeros(N,3); % point source


%% Simulation
close all;
figure;

for i = 1:time
    step = (2 * randi([0,1],N,3)) - 1; % binary (-1 or 1)
    %step = (2 * rand(N,3)) - 1;       % uniform (-1 to 1)
    %step = randn(N,3);                % std. normal (-1 to 1)
    pos = pos + L*step;
    pos(pos > bound) = bound;  
    pos(pos < -bound) = -bound;
        
    % real time particle plot
    scatter3(pos(:,1),pos(:,2),pos(:,3),'.');
    axis([-1.1*bound,1.1*bound,-1.1*bound,1.1*bound,-1.1*bound,1.1*bound]);
    label = sprintf('t = %d s',i);
    title(label);
    daspect([1 1 1]);
    drawnow;
    
end

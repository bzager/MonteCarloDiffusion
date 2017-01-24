function [pos] = absorbing(N,L,time,dt,bound,distr,init)
% Diffusion simulation with 
% parameters:
%   N -  number of particles
%   L - step size for particle
%   time - length of simulation in seconds
%   dt - size of time-step
%   bound - side length of entire domain
%   distr - probability distribution for steps
%   init - initial positions of particles

% determine prob distr for each step
if (strcmp(distr,'bin'))      % binary distribution
    step = @() (2*randi([0,1],N,2) - 1) + drift; 
elseif (strcmp(distr,'uni'))  % uniform (-1 to 1)   
    step = @() (2*rand(N,2) - 1) + drift;  
elseif (strcmp(distr,'norm')) % std. normal
    step = @() randn(N,2) + drift; 
end

pos = init;

close all;
figure;
x = bound*[1 -1 -1 1 1];
y = bound*[1 1 -1 -1 1];

numBins = 40;

for i = 1:time/dt
    prev = pos;
    
    pos = pos + L*step();
    
    pos(pos >= bound) = bound-1;   % boundary check for outer boundary
    pos(pos <= -bound) = -bound+1;
    
    % real time particle plot
    subplot(121);
    plot(x,y,'r-','LineWidth',1);
    hold on
    p1 = plot(pos(:,1),pos(:,2),'.','MarkerSize',1,'Color','k');
    axis(bound*[-1.1,1.1,-1.1,1.1]);
    label = sprintf('Particle Position at t = %.5f s',i*dt);
    title(label,'Interpreter','latex');
    xlabel('x (nm)','Interpreter','latex');
    ylabel('y (nm)','Interpreter','latex');
    daspect([1 1 1]);
    drawnow
    hold off
    
    subplot(122);
    h = histogram2(pos(:,1),pos(:,2),numBins,'Normalization','countdensity');
    label = sprintf('Particle Concentration at t=%.5f s',i*dt);
    title(label,'Interpreter','latex');
    axis([-bound bound -bound bound 0 0.002]);
    drawnow
    %[] = histcounts2();
   
end



end

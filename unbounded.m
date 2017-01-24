function [pos] = unbounded(N,L,time,dt,distr,init,a,k)
% Diffusion simulation with square boundaries
% parameters:
%   N -  number of particles   
%   time - length of simulation in seconds
%   distr - probability distribution for random walk, 'bin', 'uni', 'norm'
%   init - initial conditions of particles

drift = [a*ones(N,1) zeros(N,1)];
restore = dt*[k(1)*ones(N,1) k(2)*ones(N,1)];

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

numBins = 40; % number of intervals for histogram

for i = 1:time/dt
    pos = pos + L*step() - restore.*pos;
    
    dist = sqrt(pos(:,1).^2 + pos(:,2).^2);
    
    % real time particle plot
    subplot(121);
    plot(pos(:,1),pos(:,2),'.','MarkerSize',1,'Color','k');
    axis(1.0*10^(5)*[-1.1,1.1,-1.1,1.1]);
    label = sprintf('Particle Position at t = %.5f s',i*dt);
    title(label,'Interpreter','latex');
    xlabel('x (nm)','Interpreter','latex');
    ylabel('y (nm)','Interpreter','latex');
    daspect([1 1 1]);
    drawnow
    
    subplot(122);
    h = histogram2(pos(:,1),pos(:,2),numBins,'Normalization','countdensity');
    axis([-1.5*10^(5),1.5*10^(5),-1.5*10^(5),1.5*10^(5),0,0.1*10^-4]);
    label = sprintf('Particle Concentration at t=%.5f s',i*dt);
    title(label,'Interpreter','latex');
    drawnow
    %[] = histcounts2();
   
end

end

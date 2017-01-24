function [pos,V] = channel(N,L,time,dt,bound,wid,len,distr,init,a,U)
% Monte Carlo simulation of ion channel diffusion with membrane potential
% Compare results to Nernst relation
% parameters:
%   N -  number of particles
%   L - step size for particle
%   time - length of simulation in seconds
%   dt - size of time-step
%   bound - side length of entire domain
%   wid - width of the ion channel
%   len - length of the ion channel
%   distr - probability distribution for steps
%   init - initial positions of particles
%   a - membrane potential
%   U - resting potential
% 

% determine prob distr for each step
if (strcmp(distr,'bin'))      % binary distribution
    step = @() 2*randi([0,1],N,2) - 1; 
elseif (strcmp(distr,'uni'))  % uniform (-1 to 1)   
    step = @() 2*rand(N,2) - 1;
elseif (strcmp(distr,'norm')) % std. normal
    step = @() randn(N,2); 
end

pos = init;
a0 = a; % store initial potential
numBins = 40;
V = zeros(1,time/dt);

close all;
%figure;

% points for drawing  outer boundary
x = bound*[1 -1 -1 1 1];
y = bound*[1 1 -1 -1 1];

% points defining upper exclusion region
xu = 0.5*len*bound*[-1 1 1 -1 -1];
yu = bound*[1 1 wid/2 wid/2 1];

% points defining lower exclusion region
xl = xu;
yl = -bound*[wid/2 wid/2 1 1 wid/2];

for i = 1:time/dt
    drift = [a*ones(N,1) zeros(N,1)]; % set drift force
    prev = pos; % save last positions
    
    pos = pos + L*(step() - drift);
    
    pos(pos(:,1) <= -bound,1) = -bound ; % left boundary
    pos(pos(:,1) >= bound,1) = bound; % right boundary
    pos(pos(:,2) <= -bound,2) = -bound; % top boundary 
    pos(pos(:,2) >= bound,2) = bound; % bottom boundary
    
    outL = inpolygon(pos(:,1),pos(:,2),xl,yl); % lower exclusion region check
    outU = inpolygon(pos(:,1),pos(:,2),xu,yu); % upper exclusion region check
    
    pos(outL,:) = prev(outL,:); % return to prev pos if over bounds
    pos(outU,:) = prev(outU,:); % 
    
    
    left = sum(pos(:,1) < -len/2); % number of particles left of membrane
    right = sum(pos(:,1) > len/2); % number of particles right of membrane
     
    a = a0*(right-left)/N; % drift varies with concentration gradient
    
    V(i) = a(1);
    
    % debugging
    if (mod(i,10) == 0)
        fprintf('time = %g\n',i*dt);
        fprintf('a = %.4f\n',a);
        fprintf('right = %d\n',right);
        fprintf('left = %d\n\n',left);
        if(abs(right - left)/N < 0.01)
            fprintf('EQUILIBRIUM');
            disp(i*dt);
            return
        end
    end
    
    
    %real time particle plotting
    subplot(121);
    plot(x,y,'r-','LineWidth',1)
    patch(xu,yu,'y');
    patch(xl,yl,'y')
    hold on
    plot(pos(:,1),pos(:,2),'.','MarkerSize',1,'Color','k');
    axis(bound*[-1.1,1.1,-1.1,1.1]);
    label = sprintf('Particle Position at t = %.5f s',i*dt);
    title(label,'Interpreter','latex');
    xlabel('x (nm)','Interpreter','latex');
    ylabel('y (nm)','Interpreter','latex');
    daspect([1 1 1]);
    drawnow
    hold off
%     
    subplot(122);
    h = histogram2(pos(:,1),pos(:,2),numBins,'Normalization','countdensity','FaceAlpha',0.8);
    label = sprintf('Particle Concentration at t=%.5f s',i*dt);
    title(label,'Interpreter','latex');
    axis([-bound bound -bound bound 0 0.0002]);
    drawnow
    %[] = histcounts2();
    
    
end
end


%% parameters, grid set up, initial conditions

restart = 1;  %toggle to 1 if you want to restart the simulation, toggle to 0 if you want to continue the simulation from where it ended previously
lhorizontal = 0; %toggle to 1 for right drift
rhorizontal = 1; %toggle to 1 for left drift
vertical = 1; %toggle to 1 for vertical drift
numtimes=   6; %changing length of simulation
tmax =120*numtimes;

k=.1; %growth rate 
gamma = .0005; %diffusion coefficient
q=1/30; %drift coeficcient

n = 50;  % number of colors
beige    = [0.96 0.93 0.85];
marigold = [1.00 0.65 0.00];  % golden-orange yellow
maroon   = [0.50 0.00 0.00];

colors = [beige;
          marigold;
          maroon];

% Interpolate colormap
x  = linspace(0,1,size(colors,1));
xi = linspace(0,1,n);
cmap = interp1(x, colors, xi);


if restart == 1
    N = 18*3; %grid size
    x = [-1:2/N:1];
    dx = 2/N; %spatial discretization
    [xx,yy] = meshgrid(x,x); %create grid
    dt = (1/(15*gamma))*dx^2; %time discretization
   
    env = [ones(N/6,N+1);
    zeros(N/6,N+1);
   zeros(N/6,N+1);
     zeros(N/6,N+1);
    zeros(N/6,N+1);
     zeros(N/6+1,N+1)]; %environment (broken up into 6 columns)
    env = env';
    u = [  zeros(N/6,N/3) [
       [zeros(N/18,N/3); ...
      [zeros(N/18,N/9) ones(N/18,N/9) zeros(N/18, N/9)]; ...
      zeros(N/18,N/3)] zeros(N/6,N/3 +1) ;
    ];
        zeros(N/6,N+1);
       zeros(N/6,N+1);
        zeros(N/6,N+1);
        zeros(N/6,N+1);
        zeros(N/6+1,N+1)]; %initial distribution -- population in center of column 1
    u = u';
    
    
    disturbance = exp(-(100*(x+.9).^4)).*ones(N+1,N+1); %vertical drift
    %[ones(N/6,N+1);
       %  zeros(N/6,N+1);
       % zeros(N/6,N+1);
       %  zeros(N/6,N+1);
       %  zeros(N/6,N+1);
       %  zeros(N/6+1,N+1)]';
    
%% finite differences    
    % UXX
    uxx = zeros(N+1,N+1);
    for i = 2:N
        for j = 1:N+1
            if i == j+1
                uxx(i,j) = 1/(dx^2);
            elseif i==j
                uxx(i,j) = -(2/dx^2); %centered differences for diffusion
            elseif i == j-1
                uxx(i,j) = 1/(dx^2);
            end
        end
    end
    uxx(1,1) = -1/(dx^2);
    uxx(1,2) = 1/(dx^2);
    uxx(N+1,N+1) = -1/(dx^2);
    uxx(N+1,N) = 1/(dx^2);
    
    % UX
    if lhorizontal == 1 || vertical == 1 
        ux = zeros(N+1,N+1);
        for i = 2:N
            for j = 1:N+1
                if i == j+1
                    ux(i,j) = 1/dx; %discretization for right and vertical drift
                elseif i==j
                    ux(i,j)= -1/dx;
                end
            end
        end
        ux(1,1)=-1/dx;
        ux(N+1,N)=1/dx;
    else
        ux = zeros(N+1,N+1);
        for i = 2:N
            for j = 1:N+1
                if i == j-1
                    ux(i,j) = 1/dx; %discretization of left drift
                elseif i==j
                    ux(i,j)= -1/dx;
                end
            end
        end
        ux(1,2)=1/dx;
        ux(N+1,N+1)=-1/dx;
    end
    
    
    d2x = kron(uxx,eye(N+1));
    d2y = kron(eye(N+1),uxx);
    
    dx = kron(ux,eye(N+1));
    dy = kron(eye(N+1),ux);
end


density = zeros(1,round(tmax/(10*dt)+1));
count = 0;

%% time step 
for i = 0:dt:tmax
    count = count+1;
    if vertical == 1
        %vert = q*(disturbance.*reshape(dy*u(:),N+1,N+1));
        unew = u + ...
            dt*(gamma*(reshape(d2x*(u(:).^2),N+1,N+1)+ reshape(d2y*(u(:).^2),N+1,N+1))...
             + vert(q,i, dy, N, u, disturbance,numtimes)...
            + k*u.*(env - u));
    elseif lhorizontal == 1 || rhorizontal == 1
        unew = u + ...
            dt*(gamma*(reshape(d2x*(u(:).^2),N+1,N+1)+ reshape(d2y*(u(:).^2),N+1,N+1))...
            + horiz(q,i, dx, N, u,numtimes) ...
            + k*u.*(env - u));
    else
        unew = u + ...
            dt*(gamma*(reshape(d2x*(u(:).^2),N+1,N+1)+ reshape(d2y*(u(:).^2),N+1,N+1))...
            + k*u.*(env - u));
    end
    % if mod(count,10)==0
    %     density(count/10+1) = sum(sum(u))*(2/N)^2;
    % end
    u = unew;
    figure(1)
    pcolor(unew);
    caxis([0 1])
    colormap(cmap)
    xlabel('Column')
    ylabel("Row")
    
    title(i/(tmax/120))

    %figure(2)
    %contour(unew);
    %ylim([0 1])
    drawnow

end
sum(sum(u))*(2/N)^2

%figure(2)
%plot(density)

%% functions

%vertical drift
function out = vert(q,t,dy, N,u, disturbance,numtimes)
    if mod(t,numtimes*4)<=1*numtimes && t>4*numtimes 
    %if (t >= numtimes*23 && t<numtimes*25) || (t >= numtimes*47 && t<numtimes*49) || (t >= numtimes*71 && t<numtimes*73) || (t>= numtimes*95 && t < numtimes*97)
        out = q*(disturbance.*reshape(dy*u(:),N+1,N+1));
    else 
        out = 0;
    end
end

%horizontal drift
function out = horiz(q,t, dx, N, u,numtimes)
    %if t>0
    %if (t >= numtimes*24 && t<numtimes*25) || (t >= numtimes*48 && t< numtimes*49) || (t >= numtimes*72 && t<numtimes*73) || (t>= numtimes*96 && t < numtimes*97)
    if mod(t,numtimes*4)<=1*numtimes  && t > 4*numtimes && t < 120*numtimes-24
    %if (t >= numtimes*23 && t<numtimes*25) || (t >= numtimes*47 && t<numtimes*49) || (t >= numtimes*71 && t<numtimes*73) || (t>= numtimes*95 && t < numtimes*97) 
        out = q*(reshape(dx*u(:),N+1,N+1));
    else
        out = 0;
    end
end
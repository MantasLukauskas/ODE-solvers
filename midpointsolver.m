%% Midpoint Solver
function [t,data] = midpointsolver(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
data(1,:) = y;
for i =1:Nsteps
    dy = feval(derivs_Handle,time,y); %% evaluate the initial derivatives
    yH = y + dy*dt/2; %% take Euler step to midpoint
    dy = feval(derivs_Handle,time,y); %% re-evaluate the derivs
    y = y + dy*dt; %% shoot across the interval
    time = time+dt; %% increment time
    t(i+1) = time; %% store for output
    data(i+1,:) = y; %% store for output
end
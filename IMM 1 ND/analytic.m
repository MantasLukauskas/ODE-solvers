%% Euler Solver
%% Place this code in a file called eulersolver.m
function [t,data] = analytic(y,dt,t_final)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take.
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
data(1,:) = y;
for i =1:Nsteps
    time = time+dt;
    t(i+1) = time;
    data(i+1,:) = 0.05*(exp(-4*time)-exp(-16*time));
end
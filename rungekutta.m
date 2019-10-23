function [t,data] = rungekutta(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take.
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
data(1,:) = y;
fprintf('Step 0: t = %6.4f, w = %18.15f\n', t, y);
for i=1:Nsteps
    k1 = dt*feval(derivs_Handle,time,y);
    k2 = dt*feval(derivs_Handle,time+dt/2, y+k1/2);
    k3 = dt*feval(derivs_Handle,time+dt/2, y+k2/2);
    k4 = dt*feval(derivs_Handle,time+dt, y+k3);
    y = y + (k1+2*k2+2*k3+k4)/6;
    time = time + dt;
    t(i+1) = time;
    data(i+1,:) = y;
    fprintf('Step %d: t = %6.4f, w = %18.15f\n', i, time, y);
end
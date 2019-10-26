function [t,data] = adams3rd(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take.
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
data(1,:) = y;
fprintf('Step 0: t = %6.4f, w = %18.15f\n', t, y);

F = @(t) 0.05*(exp(-4*t)-exp(-16*t)); % true solution 
F2 = @(t) 4*0.05*exp(-4*t)-16*0.05*exp(-16*t); % true solution 

y_n = y;
y_n_1 = [F2(t(2)),F(t(2))];
y_n_2 = [F2(t(3)),F(t(3))];

for i=1:Nsteps
    k1 = 23*feval(derivs_Handle,time,y_n_2);
    k2 = 16*feval(derivs_Handle,time,y_n_1);
    k3 = 5*feval(derivs_Handle,time,y_n);
    y = y_n + dt*(k1-k2+k3)/12;
    
    time = time + dt;
    t(i+1) = time;
    data(i+1,:) = y;
    
    y_n_2 = y_n_1;
    y_n_1 = y_n;
    y_n = y;
    
    fprintf('Step %d: t = %6.4f, w = %18.15f\n', i, time, y);
end
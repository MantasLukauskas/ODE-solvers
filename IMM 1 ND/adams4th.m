function [t,data] = adams4th(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take.
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
data(1,:) = y;
fprintf('Step 0: t = %6.4f, w = %18.15f\n', t, y);

F = @(t) 0.05*(exp(-4*t)-exp(-16*t)); % true solution 
F2 = @(t) -4*0.05*exp(-4*t)+16*0.05*exp(-16*t); % true solution 

t(2) = dt
t(3) = 2*dt
t(4) = 3*dt

data(1,:) = y;
data(2,:) = [F(t(2));F2(t(2))];
data(3,:) = [F(t(3));F2(t(3))];

y_n_3 = [F(t(4));F2(t(4))];
data(4,:) = y_n_3;

for i=4:Nsteps
    
    k1 = 55*feval(derivs_Handle,t(i),data(i,:));
    k2 = 59*feval(derivs_Handle,t(i-1),data(i-1,:));
    k3 = 37*feval(derivs_Handle,t(i-2),data(i-2,:));
    k4 = 9*feval(derivs_Handle,t(i-3),data(i-3,:));
    
    y = y_n_3 + dt/24*(k1-k2+k3-k4);
    
    t(i+1) = i*dt;
    data(i+1,:) = y;
    
    y_n_3 = y;
    
    fprintf('Step %d: t = %6.4f, w = %18.15f\n', i, time, y);
end
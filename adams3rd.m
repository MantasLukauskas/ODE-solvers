function [t,data] = adams3rd(y,dt,t_final,derivs_Handle)
time = 0;
Nsteps = round(t_final/dt); %% number of steps to take.
t = zeros(Nsteps,1);
data = zeros(Nsteps,length(y));
t(1) = time; %% store intial condition
fprintf('Step 0: t = %6.4f, w = %18.15f\n', t, y);

F = @(t) 0.05*(exp(-4*t)-exp(-16*t)); % true solution 
F2 = @(t) -4*0.05*exp(-4*t)+16*0.05*exp(-16*t); % true solution 

t(2) = dt
t(3) = 2*dt

data(1,:) = y;

data(2,:) = [F(t(2));F2(t(2))];

y_n_2 = [F(t(3));F2(t(3))];
data(3,:) = [F(t(3));F2(t(3))];

fprintf("Y_2 N %d \n", y_n_2);

for i=3:Nsteps
    k1 = 23*feval(derivs_Handle,t(i),data(i,:));
    k2 = 16*feval(derivs_Handle,t(i-1),data(i-1,:));
    k3 = 5*feval(derivs_Handle,t(i-2),data(i-2,:));
    y = y_n_2 + dt/12*(k1-k2+k3);
    
    t(i+1) = i*dt;
    data(i+1,:) = y;
    y_n_2 = y;
    
    fprintf('Step %d: t = %6.4f, w = %18.15f\n', i, t(i+1), y);
end
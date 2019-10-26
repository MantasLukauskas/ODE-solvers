% To solve y' = t^-2 (sin(2t) - 2ty) s.t. y(1) = 2 for 1 <= t <= 2 using % the 3-step Adams-Bashforth. 
function [t,data] = adam_bashforth(y,dt,t_final,derivs_Handle)

time = 0;
Nsteps = round(t_final/dt);

t = zeros(Nsteps+1,1);
t(1) = time; %% store intial condition
t(2) = dt; 

data = zeros(Nsteps+1,length(y));

F = @(t) 0.05*(exp(-4*t)-exp(-16*t)); % true solution 

data(1,:) = y;

y_minus = y;
y_current = [0,F(t(2))];

% Main loop for marching N steps: 
for i = 2:Nsteps
  t(i+1) = i*dt; % time points
  
  k1 = 23*feval(derivs_Handle,t(i), [0,0]);
  k2 = 16*feval(derivs_Handle,t(i-1), y);
  %y(i+1) = y(i) + (dt/12)*(k1 - k2); %+ 5*f(t(i-2), y(i-2)) 
  
  y_further = y_current + (dt/12)*(k1 - k2);
  
  y_minus = y_current;
  y_current = y_further;

  data(i+1,:) = y_further;
  
  fprintf("Y value %d \n", y);
end 
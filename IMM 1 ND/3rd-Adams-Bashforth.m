% To solve y' = t^-2 (sin(2t) - 2ty) s.t. y(1) = 2 for 1 <= t <= 2 using % the 3-step Adams-Bashforth. 

clear, clc, clf 
f = @(t, y) (-40*y-128*y)/2;
h = 0.01; % time step size 
T = 1; % This is length of the time interval for which you're solving for, i.e. 2-(-1) = 3. 
N = T/h; % total number of times steps
F = @(t) 0.05*(exp(-4*t)-exp(-16*t)); % true solution 

  % Preallocations: 
t = zeros(N+1, 1); 
y = zeros(N+1, 1); 

% Initializations: 
t(1) = 1; % Initial time t_0 = 1 + 0*h. 
y(1) = 2; % Initial value of the ODE IVP 

% Compute the solution at t = 1+h by using the exact solution: 
t(2) = 1 + h; 
y(2) = F(t(2)); 

% Main loop for marching N steps: 
for i = 2:N 
  t(i+1) = 1 + i*h; % time points
    %Adams-Bashforth:
   y(i+1) = y(i) + (h/12)*(23*f(t(i), y(i)) - 16*f(t(i-1), y(i-1))); %+ 5*f(t(i-2), y(i-2)) 
end 

%Plotting: 
tt = linspace(1, 2, 1000);% sampling points for exact solution 

for i=1:length(tt) 
ff(i)=F(tt(i)); 
end 

plot(t,y,tt,ff,'k';) % plot the exact solution using the sampling values 
legend('3-step Adams-Bashforth', 'exact';) % adding legends
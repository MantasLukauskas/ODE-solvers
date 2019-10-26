%% derivative function
%% place in file derivs.m
function dy = derivs(time,y)
dy = zeros(2,1); %% initialize dy array and orient as column
dy(2) = (-128*y(1))/2; %% dv/dt = -x
dy(1) = y(2); %% dx/dt = v

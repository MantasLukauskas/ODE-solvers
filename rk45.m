function [t,data] = rk45(y,dt,t_final,epsilon,derivs_Handle)
    i = 0;
    q = 0;
    time = 0;
    t(1) = time; %% store intial condition
    data(1,:) = y;
    fprintf("Step %d: t = %6.4f, w = %18.15f\n", i, time, y);
    while time < t_final
        q = q + 1;
        %fprintf("Step %d\n", q);
        dt = min(dt, 2-time);
        k1 = dt*feval(derivs_Handle,time,y);
        k2 = dt*feval(derivs_Handle,time+dt/4, y+k1/4);
        k3 = dt*feval(derivs_Handle,time+3*dt/8, y+3*k1/32 + 9*k2/32);
        k4 = dt*feval(derivs_Handle,time+12*dt/13, y+1932*k1/2197-7200*k2/2197+7286*k3/2197);
        k5 = dt*feval(derivs_Handle,time+dt, y+439*k1/216 - 8*k2+3680*k3/513-845*k4/4104);
        k6 = dt*feval(derivs_Handle,time+dt/2, y - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40);
        w1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
        w2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55;
        
        R = abs(w1-w2)/dt;
                
        delta = 0.84*(epsilon/R(2))^(1/4);
        
        if R <= epsilon
            time = time+dt;
            t(i+1) = time;
            y = w1;
            data(i+1,:) = y;
            i = i+1;
            fprintf("Step %d: t = %6.4f, w = %18.15f\n", i, time, y);
            dt = delta*dt;
        else
            dt = delta*dt;
        end
    end
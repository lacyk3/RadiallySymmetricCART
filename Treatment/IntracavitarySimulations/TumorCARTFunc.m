function [FG] = TumorCARTFunc(u_n, u_nn, dt, dr, DT, DC, beta_T, gamma, s,l, alpha, zeta, beta_C, chi, r, u_star)

v_nn = u_nn(length(u_nn)/2 + 1:end);
v_n = u_n(length(u_n)/2 + 1:end);
u_nn = u_nn(1:length(u_nn)/2);
u_n = u_n(1:length(u_n)/2);

F = 0*u_n;
G = 0*v_n;
F1_u_nn = F;
F1_u_n = F;
F2_v_nn = G;
F2_v_n = G;

if u_n(1) > 0 && v_n(1) > 0
F1_u_nn(1) = (1 - beta_T*u_nn(1) - gamma*v_nn(1)^l/(s*u_nn(1)^l + v_nn(1)^l))*u_nn(1);
F1_u_n(1) = (1 - beta_T*u_n(1) - gamma*v_n(1)^l/(s*u_n(1)^l + v_n(1)^l))*u_n(1);

F2_v_nn(1) = alpha*(v_nn(1)^(2*l)*u_nn(1)^2/((s*u_nn(1)^l + v_nn(1)^l)^2 + zeta*v_nn(1)^(2*l)*u_nn(1)^2) - beta_C*u_nn(1) - chi)*v_nn(1);
F2_v_n(1) = alpha*(v_n(1)^(2*l)*u_n(1)^2/((s*u_n(1)^l + v_n(1)^l)^2 + zeta*v_n(1)^(2*l)*u_n(1)^2) - beta_C*u_n(1) - chi)*v_n(1);
end

if u_n(1) <= 0 && v_n(1) > 0
F2_v_nn(1) = -alpha*chi*v_nn(1);
F2_v_n(1) = -alpha*chi*v_n(1);
end


F(1) = u_nn(1) - (dt/2)*(DT*(2*u_nn(2) - 2*u_nn(1))/dr^2 + F1_u_nn(1)) - u_n(1) - (dt/2)*(DT*(2*u_n(2) - 2*u_n(1))/dr^2 + F1_u_n(1));
G(1) = v_nn(1) - (dt/2)*(DC*(2*v_nn(2) - 2*v_nn(1))/dr^2 + F2_v_nn(1)) - v_n(1) - (dt/2)*(DC*(2*v_n(2) - 2*v_n(1))/dr^2 + F2_v_n(1));


if u_n(end) > 0 && v_n(end) > 0
F1_u_nn(end) = (1 - beta_T*u_nn(end) - gamma*v_nn(end)^l/(s*u_nn(end)^l + v_nn(end)^l))*u_nn(end);
F1_u_n(end) = (1 - beta_T*u_n(end) - gamma*v_n(end)^l/(s*u_n(end)^l + v_n(end)^l))*u_n(end);

F2_v_nn(end) = alpha*(v_nn(end)^(2*l)*u_nn(end)^2/((s*u_nn(end)^l + v_nn(end)^l)^2 + zeta*v_nn(end)^(2*l)*u_nn(end)^2) - beta_C*u_nn(end) - chi)*v_nn(end);
F2_v_n(end) = alpha*(v_n(end)^(2*l)*u_n(end)^2/((s*u_n(end)^l + v_n(end)^l)^2 + zeta*v_n(end)^(2*l)*u_n(end)^2) - beta_C*u_n(end) - chi)*v_n(end);
end

if u_n(end) <= 0 && v_n(end) > 0
F2_v_nn(end) = -alpha*chi*v_nn(end);
F2_v_n(end) = -alpha*chi*v_n(end);
end

F(end) = u_nn(end) - (dt/2)*(DT*(2*u_nn(end-1) - 2*u_nn(end))/dr^2 + F1_u_nn(end)) - u_n(end) - (dt/2)*(DT*(2*u_n(end-1) - 2*u_n(end))/dr^2 + F1_u_n(end));
G(end) = v_nn(end) - (dt/2)*(DC*(2*v_nn(end-1) - 2*v_nn(end))/dr^2 + F2_v_nn(end)) - v_n(end) - (dt/2)*(DC*(2*v_n(end-1) - 2*v_n(end))/dr^2 + F2_v_n(end));

for m = 2:length(u_n)-1
    
    if u_n(m) > 0 && v_n(m) > 0
        F1_u_nn(m) = (1 - beta_T*u_nn(m) - gamma*v_nn(m)^l/(s*u_nn(m)^l + v_nn(m)^l))*u_nn(m);
        F1_u_n(m) = (1 - beta_T*u_n(m) - gamma*v_n(m)^l/(s*u_n(m)^l + v_n(m)^l))*u_n(m);
 
        F2_v_nn(m) = alpha*(v_nn(m)^(2*l)*u_nn(m)^2/((s*u_nn(m)^l + v_nn(m)^l)^2 + zeta*v_nn(m)^(2*l)*u_nn(m)^2) - beta_C*u_nn(m) - chi)*v_nn(m);
        F2_v_n(m) = alpha*(v_n(m)^(2*l)*u_n(m)^2/((s*u_n(m)^l + v_n(m)^l)^2 + zeta*v_n(m)^(2*l)*u_n(m)^2) - beta_C*u_n(m) - chi)*v_n(m);
    end
    
    if u_n(m) <= 0 && v_n(m) > 0
        F2_v_nn(m) = -alpha*chi*v_nn(m);
        F2_v_n(m) = -alpha*chi*v_n(m);
    end
    
        if u_n(m+1) > u_star && u_n(m-1) > u_star
            F(m) = u_nn(m) - u_n(m) - (dt*DT/r(m))*(u_nn(m+1) - u_nn(m-1) + u_n(m+1) - u_n(m-1))/(2*dr) - (dt/2)*(F1_u_n(m) + F1_u_nn(m)) - ...
                   dt*DT*(u_nn(m+1) - 2*u_nn(m) + u_nn(m-1) + u_n(m+1) - 2*u_n(m) + u_n(m-1))/(2*dr^2);
        elseif u_n(m+1) <= u_star && u_n(m-1) > u_star
            F(m) = u_nn(m) - (dt/2)*(DT*(2*u_nn(m-1) - 2*u_nn(m))/dr^2 + F1_u_nn(m)) - u_n(m) - (dt/2)*(DT*(2*u_n(m-1) - 2*u_n(m))/dr^2 + F1_u_n(m));
        elseif u_n(m-1) <= u_star && u_n(m+1) > u_star
            F(m) = u_nn(m) - (dt/2)*(DT*(2*u_nn(m+1) - 2*u_nn(m))/dr^2 + F1_u_nn(m)) - u_n(m) - (dt/2)*(DT*(2*u_n(m+1) - 2*u_n(m))/dr^2 + F1_u_n(m));
        else
            F(m) = u_nn(m) - u_n(m) - (dt/2)*(F1_u_n(m) + F1_u_nn(m));
        end
        
        G(m) = v_nn(m) - v_n(m) - (dt*DC/r(m))*(v_nn(m+1) - v_nn(m-1) + v_n(m+1) - v_n(m-1))/(2*dr) - (dt/2)*(F2_v_n(m) + F2_v_nn(m)) - ...
                   dt*DC*(v_nn(m+1) - 2*v_nn(m) + v_nn(m-1) + v_n(m+1) - 2*v_n(m) + v_n(m-1))/(2*dr^2);
        
end

FG = [F; G];

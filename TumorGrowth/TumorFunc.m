function [F] = TumorCARTFunc(u_n, u_nn, dt, dr, DT, beta_T, r, u_star)


F = 0*u_n;

F1_u_nn = F;
F1_u_n = F;



F1_u_nn(1) = (1 - beta_T*u_nn(1))*u_nn(1);
F1_u_n(1) = (1 - beta_T*u_n(1))*u_n(1);




F(1) = u_nn(1) - (dt/2)*(DT*(2*u_nn(2) - 2*u_nn(1))/dr^2 + F1_u_nn(1)) - u_n(1) - (dt/2)*(DT*(2*u_n(2) - 2*u_n(1))/dr^2 + F1_u_n(1));

F1_u_nn(end) = (1 - beta_T*u_nn(end))*u_nn(end);
F1_u_n(end) = (1 - beta_T*u_n(end))*u_n(end);




F(end) = u_nn(end) - (dt/2)*(DT*(2*u_nn(end-1) - 2*u_nn(end))/dr^2 + F1_u_nn(end)) - u_n(end) - (dt/2)*(DT*(2*u_n(end-1) - 2*u_n(end))/dr^2 + F1_u_n(end));

for m = 2:length(u_n)-1
    
    
        F1_u_nn(m) = (1 - beta_T*u_nn(m))*u_nn(m);
        F1_u_n(m) = (1 - beta_T*u_n(m))*u_n(m);


    
    
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
        
       
        
end
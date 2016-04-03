function [T_stored] = Store_temp_loss(T_stored, r_store, L_store, ...
    lambda_ins, ins_th, st_time, n_stored, c_air_molar, T_am)

R_th_cy = log((r_store + ins_th)/r_store)*(1/(2*pi*L_store*lambda_ins));
A_end = pi * (r_store^2);
R_th_end = (ins_th)/(A_end*lambda_ins);

dt = 10; time_index = st_time*3600/dt;

for i = 1:time_index
    
rate_loss = ((T_stored-T_am)/R_th_cy)+((T_stored-T_am)/R_th_end);
loss = dt*rate_loss;

T_stored = (n_stored*c_air_molar*T_stored - loss)/(n_stored*c_air_molar);

if T_stored<T_am
    T_stored = T_am;
end

end



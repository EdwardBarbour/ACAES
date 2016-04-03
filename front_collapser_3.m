function [T_peb, Heat_Ex_loss, T_PB_anim_test] = front_collapser_3(T_peb, k_peb, L_bed, r_bed, c_peb, dl, ...
    st_time, mass_per_dh, ins_th, lambda_ins, T_am)

dh = dl;
dt = 10;
A = (r_bed^2)*pi;
time_index = (st_time*3600)/dt;
PB_index = L_bed/dl;

%thermal resistance of slice of the cylinder
R_th_cy = log((r_bed + ins_th)/r_bed)*(1/(2*pi*dl*lambda_ins));

%thermal resistance of end
R_th_end = (ins_th)/(A*lambda_ins);
Heat_Ex_loss = 0;

%%%%% for animation %%%%%%%
interval = 300/dt; %% every 5 min save temp profile
T_PB_anim_test = zeros(PB_index, time_index/interval);
x = 1;

for j = 1:time_index;
for i = 1:PB_index
if(i==1)
    rate_loss = ((T_peb(i)-T_am)/R_th_cy)+((T_peb(i)-T_am)/R_th_end);
    Loss = rate_loss*dt;
    T_peb(i) = (mass_per_dh*c_peb*T_peb(i) - Loss)/(mass_per_dh*c_peb); 
    if(T_peb(i)<T_am);T_peb(i)=T_am;end
else
        if(i==PB_index)
            rate_heat_gain = (((T_peb(i-1)-T_peb(i))/dh)*k_peb*A);
            rate_loss = ((T_peb(i)-T_am)/R_th_cy)+((T_peb(i)-T_am)/R_th_end);
            Loss = rate_loss*dt;
            heat_gain = rate_heat_gain*dt;
        else
            rate_heat_gain = (((T_peb(i-1)-T_peb(i))/dh)*k_peb*A);
            rate_loss = ((T_peb(i)-T_am)/R_th_cy);
            Loss = rate_loss*dt;
            heat_gain = rate_heat_gain*dt;
        end
        
    T_peb_1 = T_peb(i-1); T_peb_2 = T_peb(i);
    if(Loss>mass_per_dh*c_peb*(T_peb(i)-T_am));Loss=mass_per_dh*c_peb*(T_peb(i)-T_am);end    
    Q_peb_1 = mass_per_dh*c_peb*T_peb(i-1); Q_peb_2 = mass_per_dh*c_peb*T_peb(i);  
    T_peb(i-1) = (mass_per_dh*c_peb*T_peb(i-1) - heat_gain)/(mass_per_dh*c_peb);
    T_peb(i) = (mass_per_dh*c_peb*T_peb(i) + heat_gain - Loss)/(mass_per_dh*c_peb);
    
    T_limit = (Q_peb_1+Q_peb_2-Loss)/(2*c_peb*mass_per_dh);

    if(T_peb_1>T_peb_2)
    if(T_peb(i)>T_limit)
        T_peb(i)=T_limit; T_peb(i-1) = T_limit; 
    end
    end

    if(T_peb_1<=T_peb_2)
    if(T_peb(i)<T_limit)
        T_peb(i)=T_limit; T_peb(i-1)=T_limit;
    end
    end
    
    if(T_peb(i)<T_am);T_peb(i)=T_am;end 
    if(T_peb(i-1)<T_am);T_peb(i-1)=T_am;end
    
    
    
% test=isnan(T_peb(i-1));
test = T_peb(i-1);
if(test==inf)
fprintf('%f \n', test)
fprintf('%f \n', T_peb(i-1))
fprintf('%f \n', T_peb(i-2))
fprintf('%f \n', T_peb(i-3))
fprintf('%f \n', heat_gain)
fprintf('%f \n', Loss)
error('myApp:argChk', 'nan')
end

end

% Heat_Ex_loss = Heat_Ex_loss + Loss*(1 - (T_am/T_peb(i)));
Heat_Ex_loss = Heat_Ex_loss + c_peb*mass_per_dh*(Loss/(c_peb*mass_per_dh) - T_am*log((T_peb(i)+(Loss/(c_peb*mass_per_dh)))/T_peb(i)));
    
end

% return temperature profile for animation;

if(rem(j,interval)==0)
T_PB_anim_test(:,x) = T_peb;
x = x + 1;
end

end


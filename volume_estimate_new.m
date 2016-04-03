function [store_volume, total_moles] = volume_estimate_new(p, T_test, store_volume, p_i, R, ...
    NS, work_required, k, eta_pol_c, c_p_air)

T_am = T_test(1);
V_s = store_volume;
Work = zeros(NS, 1);

for stage = 1:NS
    if(p(2*stage+1)<p_i)
        W = (c_p_air*(p(2*NS+1)-p_i)*V_s/R) * ((p(2*stage)/p(2*stage-1))^((k-1)/(k*eta_pol_c)) - 1);
    else
        if stage == NS
        if p_i>p(2*stage-1)
        W = (c_p_air*V_s*p(2*stage)/R) * ((p_i/p(2*stage)) - 1 + ((eta_pol_c*k)/(k-1+eta_pol_c*k))* ...
            ((p(2*stage)/p(2*stage-1))^((k-1)/(eta_pol_c*k)) - ...
            (p_i/p(2*stage))*((p_i/p(2*stage-1))^((k-1)/(eta_pol_c*k))) ) );            
        else
        W = (c_p_air*V_s*p(2*stage)/R) * ((p(2*stage-1)/p(2*stage)) - 1 + ((eta_pol_c*k)/(k-1+eta_pol_c*k))* ...
            ((p(2*stage)/p(2*stage-1))^((k-1)/(eta_pol_c*k)) - ...
            (p(2*stage-1)/p(2*stage)) ) );
        end
        else
        if p_i>p(2*stage-1)
        W = (c_p_air*V_s*p(2*stage)/R) * ((p_i/p(2*stage)) - 1 + ((eta_pol_c*k)/(k-1+eta_pol_c*k))* ...
            ((p(2*stage)/p(2*stage-1))^((k-1)/(eta_pol_c*k)) - ...
            (p_i/p(2*stage))*((p_i/p(2*stage-1))^((k-1)/(eta_pol_c*k)))) ) ...
            + (c_p_air*(p(2*NS+1)-p(2*stage+1))*V_s/R) * ((p(2*stage)/p(2*stage-1))^((k-1)/(k*eta_pol_c)) - 1);
        else
        W = (c_p_air*V_s*p(2*stage)/R) * ((p(2*stage-1)/p(2*stage)) - 1 + ((eta_pol_c*k)/(k-1+eta_pol_c*k))* ...
            ((p(2*stage)/p(2*stage-1))^((k-1)/(eta_pol_c*k)) - ...
            (p(2*stage-1)/p(2*stage))) ) ...
            + (c_p_air*(p(2*NS+1)-p(2*stage+1))*V_s/R) * ((p(2*stage)/p(2*stage-1))^((k-1)/(k*eta_pol_c)) - 1);    
        end
        end
    end
    Work(stage)=100000*W;
end

Total_work = sum(Work)/3600;

store_volume = store_volume*work_required/Total_work;
total_moles = (p(2*NS+1)*100000*store_volume)/(R*T_am);

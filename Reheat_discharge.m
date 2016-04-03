function[T_air_out, T_peb, T_fluid, dn_out, mass_per_dh, p_air_out, h_vol_lof, Heat_Ex_loss, p_bed, moles_in_slice, Delta_T_fluid] = ...
    Reheat_discharge(dn, N_dot, T_air_in, p_air_in, c_air_mol, T_peb, T_fluid, length_PB, r_bed, dl, c_peb, k_peb, ...
    lambda_ins, ins_th, p_store, V_store, T_store, T_am, p_bed, moles_in_slice, Delta_T_fluid, voidage)

%Specify the geometry of the PBHE
% stuff to be read in
height = length_PB;
diameter = 2*r_bed;
dh = dl;

% stuff to be defined
pebble_size = 0.01;
rho_peb = 2640;
mu = 0.00003; %air viscosity

n_slices = round(height/dh);
A = pi*(diameter/2)^2;
void_vol_slice = dh*A*voidage;
mass_pebbles = (1-voidage)*rho_peb*height*A;
mass_per_dh = mass_pebbles/n_slices;

%Lof and Hawley correlation for volumetric h
G = (N_dot*0.029)/(A);
% h_vol_lof = 652*(G/pebble_size)^0.7;
h_vol_lof = 700*(G/pebble_size)^0.76;
h_vol = h_vol_lof;

%losses to environment
%thermal resistance of slice of the cylinder
R_th_cy = log((r_bed + ins_th)/r_bed)*(1/(2*pi*dl*lambda_ins));
%thermal resistance of end
R_th_end = (ins_th)/(A*lambda_ins);

% problems/assumptions
% 1) the pressure throughout PB is always same as incoming air
% 2) the moles in each slice remains constant despite temperature changes

%assume that air mixes instantaneously as it enters slice and the dn added
%pushes a dn into the next slice

dns = zeros(1,n_slices+1);
% Calculate the dn's, as they are calculated wrt previous temps, need to be
% small compared to moles_in_slice
% if p_air > p_store then need to consider pressure increases in PB

% in order to try and ensure that the pressure in the PB remains constant
% the moles pushed out is calculated using the ESTIMATED next air
% temperature.

% moles in each slice depends on bed pressure and temp

% store the old fluid temperature to calculate the dynamic response of the
% bed
T_fluid_old = T_fluid;

if(p_air_in>=p_store)
    
    % this is different during charge compared to discharge.
    % think this now works
    dns(1) = dn;
    dns(2) = dns(1) + moles_in_slice(1)*(1-T_fluid(1)/(T_fluid(1)+Delta_T_fluid(1))) + dn*(T_store*void_vol_slice/(V_store*(T_fluid(1)+Delta_T_fluid(1))));
    for i = 3:n_slices+1
        dns(i) = dns(i-1) + moles_in_slice(i-1)*(1-T_fluid(i-1)/(T_fluid(i-1)+Delta_T_fluid(i-1))) ... 
             + dn*(T_store*void_vol_slice/(V_store*(T_fluid(i-1)+Delta_T_fluid(i-1))));
    end
    p_bed = p_bed - dn*8.31*T_store/(100000*V_store);
else
    
    % Think this is working reasonably well...
    
    % pressure is constant, only need consider Temperature differences
    dns(1) = dn;
    dns(2) = dns(1) + moles_in_slice(1)*(1-T_fluid(1)/(T_fluid(1)+Delta_T_fluid(1)));
    for i = 3:n_slices+1
        dns(i) = dns(i-1) + moles_in_slice(i-1)*(1-T_fluid(i-1)/(T_fluid(i-1)+Delta_T_fluid(i-1)));
    end

end

% Update the moles in each slice as dn(slice+1) is pushed out
moles_in_slice = moles_in_slice-dns(2:n_slices+1);

% Calculate the new fluid temperature after dn(slice) has been added
T_fluid(2:n_slices) = (moles_in_slice(2:n_slices).*T_fluid(2:n_slices) + ...
    dns(2:n_slices).*T_fluid(1:n_slices-1))./(moles_in_slice(2:n_slices)+dns(2:n_slices));
% need to have this at the end or it distorts the above calculation
T_fluid(1) = (moles_in_slice(1)*T_fluid(1) + dns(1)*T_air_in)/(moles_in_slice(1)+dns(1));

% correct the moles in each slice of the exchanger
moles_in_slice = moles_in_slice+dns(1:n_slices);
 
% approximate time spent in first slice: moles_in_slice/molar_flowrate
dt_1 = moles_in_slice(1)/N_dot;
% dt_1 is the time to replace all the moles in 1 slice
dt = (dns(1)/moles_in_slice(1))*dt_1;
clear dt_1;

%heat movement terms
Q_fluid_to_peb = h_vol*(T_fluid - T_peb)*dt;
Q_loss = ((T_peb-T_am)/R_th_cy)*dt;
Q_loss(1) = Q_loss(1)+((T_peb(1)-T_am)/R_th_end)*dt;
Q_loss(n_slices) = Q_loss(n_slices)+((T_peb(n_slices)-T_am)/R_th_end)*dt;
Q_cond = ((k_peb*A)/dh)*diff(T_peb)*dt;
%Q cond is a vector n-1 in length

% Now the new total heats (order doesn't matter for bed terms)
Q_fluid = c_air_mol*moles_in_slice.*T_fluid; Q_peb = (c_peb*mass_per_dh*T_peb);

% calculate limits on temps and maximum Q terms
% first do limit air to solid heat transfer
T_limit_fluid_to_peb = (Q_fluid+Q_peb)./(c_peb*mass_per_dh + c_air_mol*moles_in_slice);
Q_fluid_to_peb_max = (T_limit_fluid_to_peb - T_peb)*c_peb*mass_per_dh;
% check if sign of heat transfers is different..
s = sign(Q_fluid_to_peb.*Q_fluid_to_peb_max);
Q_fluid_to_peb(abs(Q_fluid_to_peb)>abs(Q_fluid_to_peb_max))=Q_fluid_to_peb_max(abs(Q_fluid_to_peb)>abs(Q_fluid_to_peb_max));
Q_fluid_to_peb = Q_fluid_to_peb.*s;
clear s;
% conduction limit
Q_cond_max = c_peb*mass_per_dh*diff(T_peb)/2;
% don't need check if sign of heat transfers is different..
% s = sign(Q_cond.*Q_cond_max);
Q_cond(abs(Q_cond)>abs(Q_cond_max))=Q_cond_max(abs(Q_cond)>abs(Q_cond_max));
% Q_cond = Q_cond.*s;
% loss limit, this time should always have the same sign.
Q_loss_max = (T_peb - T_am)*c_peb*mass_per_dh;
Q_loss(abs(Q_loss)>abs(Q_loss_max))=Q_loss_max(abs(Q_loss)>abs(Q_loss_max));

% Now calculate the new temperatures

T_fluid = (Q_fluid-Q_fluid_to_peb)./(c_air_mol*moles_in_slice);
T_peb(1) = (Q_peb(1) + Q_fluid_to_peb(1) - Q_loss(1) + Q_cond(1))/(c_peb*mass_per_dh);
T_peb(n_slices) = (Q_peb(n_slices) + Q_fluid_to_peb(n_slices) - Q_loss(n_slices) - Q_cond(n_slices-1))/(c_peb*mass_per_dh);
T_peb(2:n_slices-1) = (Q_peb(2:n_slices-1) + Q_fluid_to_peb(2:n_slices-1) - ...
    Q_loss(2:n_slices-1) + Q_cond(2:n_slices-1) - Q_cond(1:n_slices-2))/(c_peb*mass_per_dh);

% use the temperature difference as the expected temperature difference in
% the next addition of dn moles
Delta_T_fluid = T_fluid - T_fluid_old;

%%%%%%%%%%%%%%%%%%%% ERGUN EQUATION %%%%%%%%%%%%%%%%%
vol_flow_rate = (N_dot*8.31*mean(T_fluid))/(p_bed*100000);

V_sup = vol_flow_rate/A;

V_s = V_sup;
rho = (0.029*p_bed*100000)/(8.31*mean(T_fluid));
Gr_p = (pebble_size*V_s*rho)/((1-voidage)*mu);
f = 150/Gr_p + 1.75;

Delta_p = ((f*length_PB*rho*(V_s^2))/pebble_size)*((1-voidage)/(voidage^3));

% allowing for the pressure drop to reduce pressure through the PB
p_air_out = (p_bed*100000 - Delta_p)/100000;

% Re_star = (rho*V_sup*pebble_size)/(mu*(1-voidage));
% Re_star is the same as Gr_p
% Pr = ((c_air_mol/0.029)*mu)/0.04;
% Nu = 0.023*(Re_star^0.8)*(Pr^0.33);
% Nu = 2 + alpha * Re^0.5 * Pr^0.3;
% h = Nu*pebble_size/0.04;

% temperature out is the fluid temperature in the last slice
T_air_out=T_fluid(n_slices);

% Pec = Re_star*Pr;

% moles out is the moles pushed from the last slice
dn_out = dns(n_slices);

% heat exergy loss
Heat_Ex_loss = sum(  c_peb*mass_per_dh*(Q_loss/(c_peb*mass_per_dh) - T_am*log((T_peb+(Q_loss/(c_peb*mass_per_dh)))./T_peb))  );
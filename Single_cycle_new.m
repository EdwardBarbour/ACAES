% An attempt at a dynamic model of an AA-CAES system with Packed Beds to
% store the compression heat

% Edward Barbour
% edward.r.barbour@gmail.com
% 10/11/2014

clc; clear all;

% Parameters that need specified

kWh=3600000; % for converting J to kWh
store_volume = 1; %[m^3]
T_am = 290; p_am = 1.01325; % ambient atmospheric temperature [k] and pressure [bar]
store_temp = T_am; %set inital tore temperature to ambient
charge_time = 4; %hours TIME TO CHARGE
st_time = 10; %time between end of charge and start of discharge - for calcualting thermal losses
idle_time = 6; %time between end of discharge and start of next charge - for calcualting thermal losses
MM_air = 0.02897; % molar mass of air kg
R = 8.31; heat_cap_air = 1010; %kJ per kg per K
c_air_molar = heat_cap_air*MM_air; %kJ per mol per K
eta_pol_c = 0.85; eta_pol_t = 0.85; % compressor and turbine efficiencies
k = 1.4; % cp/cv exponent

% the variable used for the pressure estimation calculation
r = 1;

p_i = 20*p_am; % initial pressure
final_pressure = 80*p_am; %Choose the final pressure
NS = 2; % Specify the number of stages
work_required = 2000; %kWh specify the work required and estimate the volume required

pressure_drop = 0.02;
% the ability to guess different pressure drops for each stage
p_drop_guess(1:NS)=pressure_drop;
% for recording pressure drops
Press_d = 0; P_drop_max = zeros(1, NS); P_drop_min = zeros(1, NS); P_drop = zeros(1, NS);
% for recording h_vol
h_vol_min = 0; h_vol_max = 0;

%Packed bed characteristics
c_peb = 1000; %specific heat capacity of pebbles [J/(kgK)]
k_peb = 4; %thermal conductivity of pebbles
L_bed = [12 12]; %length of packed beds in metres [m]
r_bed = [0.6 0.6]; %radius of packed in metres
%estimate VOLUMETRIC Heat transfer coefficient [W/(m3K)].
dl = 0.02; %length of slice [m]
% insulation characteristics of the thermal stores
ins_th = 0.2; %insulation thickness [m]
lambda_ins = 0.3; %insulation thermal conductivity [W/(m·K]
voidage = 0.36; % if chnge this also need to change in Reheat program

dn_in = 0.5; % mole increment

% calculating the intermediate pressures
p_test = zeros(1,NS*2+1); T_test = zeros(1,NS*2+1); p_bed = zeros(1, NS);

% This section yields a quick iterative calculation to reach the desired
% final storage pressure.

while(p_test(NS*2+1)<final_pressure) 

a = r; b = p_test(NS*2+1);
r = r+0.1;

T_test(1) = T_am; p_test(1) = p_am;

for i = 1:NS
j = 2*i - 1;
p_test(j+1) = r*p_test(j); %pressure increased by compression
T_test(j+1) = adiabatic_comp(T_test(j), p_test(j+1), p_test(j), k, eta_pol_c); %temperature increased by compression
T_test(j+2) = HE(T_test(j+1), 1, T_am); %temperature decreased by Heat Exchanger
p_test(j+2) = p_test(j+1) - pressure_drop; %pressure decreased by pressure drop through the heat exchanger
end

c = r; d = p_test(NS*2+1);

end
diff1 = (final_pressure-b); diff2 = (-final_pressure+d);
r = a + (0.1/(diff1+diff2))*(diff1); %r is now the pressure ratio required for the simulation
clearvars a b c d 

% Now calculate the reference pressure and temps with the correct r (pressure ratio) value
T_test(1) = T_am; p_test(1) = p_am; active_stage = 0;
for i = 1:NS
j = 2*i - 1;
p_test(j+1) = r*p_test(j);
T_test(j+1) = adiabatic_comp(T_test(j), p_test(j+1), p_test(j), k, eta_pol_c);
T_test(j+2) = HE(T_test(j+1), 1, T_am); 
p_test(j+2) = p_test(j+1) - pressure_drop;
% assign the initial pressure of the packed beds and
% find the initial active compression stages
p_bed(i) = p_test(j+1); if(p_bed(i)<p_i+pressure_drop); active_stage = active_stage+1; end 
if(p_bed(i)>p_i+pressure_drop); p_bed(i)=p_i+pressure_drop; end
end

% now that reference pressures are known can calculate volume
[store_volume, total_moles] = volume_estimate_new(p_test, T_test ...
    , store_volume, p_i, R, NS, work_required, k, eta_pol_c, c_air_molar/1000);

% assume store cylinder with length 5 times radius - 
% need the store geometry to estimate the thermal loss
r_store = (store_volume/(pi*5))^(1/3); L_store = 5*r_store;

%number of cycles to be simulated - cannot be changed in this version of
%code. Just kept it for the next 16 lines
n_cy_max = 1;

% logging temperature profiles of the packed beds
% create a variable for each PB that stores temperature profile at the start of each cycle, end of compression and just before expansion 
for z = 1:NS
    z1 = num2str(z);
    name1 = (['Temp_PB_',z1,'_initial']);
    v = genvarname(name1, who);
    eval([v ' = zeros(n_cy_max, L_bed(z)/dl);']);
    name1 = (['Temp_PB_',z1,'_post_comp']);
    v = genvarname(name1, who);
    eval([v ' = zeros(n_cy_max, L_bed(z)/dl);']);
    name1 = (['Temp_PB_',z1,'_pre_exp']);
    v = genvarname(name1, who);
    eval([v ' = zeros(n_cy_max, L_bed(z)/dl);']);
end

n_cycle = 1;

% calculate the initial moles contained
n_i = (p_i*100000*store_volume)/(R*T_am); V_store = store_volume; 

%Making a cell array to contain temperature of pebbles and the air in the
%PB regenerators. Also need the dynamic paramters for the new PB model
Bed_index = L_bed/dl; Area = pi*(r_bed.*r_bed); V_bed = L_bed.*Area*voidage; void_vol_slice = dl*Area*voidage;
A = [NS, 1];
PB_peb = cell(A);
PB_fluid = cell(A);
Delta_T_fluid = cell(A);
moles_in_slice = cell(A);
Delta_n = cell(A);

for z = 1:NS
    entry1(1:Bed_index(z))=T_am; entry2(1:Bed_index(z))=0;
    PB_peb{z}=entry1;
    PB_fluid{z}=entry1;
    Delta_T_fluid{z}=entry2;
    entry3 = (p_bed(z)*100000*void_vol_slice(z))./(8.31*PB_fluid{z});
    moles_in_slice{z} = entry3;
    entry4(1:Bed_index(z)) = (dn_in/Bed_index(z))*((V_bed(z)/V_store)/(1+V_bed(z)/V_store));
    Delta_n{z} = entry4;
    clear entry1 entry2 entry3 entry4;
end
% The moles left in the packed bed must be calculated 
dn_out = dn_in - sum(Delta_n{active_stage}); clear z;

% store the temperature profiles of the bed and the air in an appropriately
% named variable. Each row corresponds to a cycle
for z = 1:NS
    z1 = num2str(z);    
    name1 = (['Temp_PB_',z1,'_initial']);
    eval([name1 '(n_cycle,:) = PB_peb{z};']);
end

%now calculate the average molar flow rate
N_dot = (total_moles-n_i)/(charge_time*3600);

% Variables for storing run information
index_max = ceil(p_test(NS*2+1));
pressure = ones(1,index_max); moles = zeros(1,index_max); Storage_temperature(1:index_max)=T_am; moles(1) = n_i;
% The variable for incrementing the logging
temp_var = 1;

%Re-initialise the working other arrays
p_air = zeros(1,NS*2+1); T_air = zeros(1,NS*2+1); p_air_ex = zeros(1,NS*2+1); T_air_ex = zeros(1,NS*2+1);
% p_air and T_air are pressure & temperature for compression
% p_air_ex and T_air_ex are pressure & temp for expansion
%Put in the initial conditions
p_stored = p_i; p_air(1) = p_am; T_air(1) = T_am;

% Initialising terms like initial moles stored, volume and work terms
V_stored = store_volume;
n_stored = n_i; T_stored = T_am; W_stored = zeros(1,NS); W_released = zeros(1,NS);
Q_stored = n_stored*c_air_molar*T_stored;

n_test = zeros(1,NS);

% Initialising exergy loss terms
Comp_Ex_dest=zeros(1,NS); Turb_Ex_dest=zeros(1,NS); Exit_loss = 0;
PB_Ex_loss_comp=zeros(1,NS); PB_Ex_loss_st=zeros(1,NS); PB_Ex_loss_turb=zeros(1,NS);

% matricies for logging pressure and temperature against store pressure and
% stages
pressure_mat_in = zeros(temp_var, NS*2+1);
temperature_mat_in = zeros(temp_var, NS*2+1);

%firstly need to specify roughly the first dn that enters the store
dn_into_store = dn_in - sum(Delta_n{active_stage});

for i = 1:NS

while(p_stored<p_test(2*i+1))
% when the pressure reaches this threshold, then i is incremented and dn
% will pass through stages 1 to i. The threshold is then increased. Initially dn only passes through i stages. 
for stage = 1:i
    
j = 2*stage - 1;
% add dn moles of air to the store
% through the compressor

if(stage == 1)
    dn = dn_in;
end

if(stage==i)  
p_stored_high = ((n_stored+dn_into_store)*R*T_stored/store_volume)/100000;
p_air(j+1) = p_stored_high + p_drop_guess(stage);
% The air is compressed to just above the pressure in the store to allow
% the air to flow in. The pressure drop term is added to make sure the air
% has enough pressure to flow in despite pressure losses in the stage before it enters the HP air store 
if(p_air(j+1)<p_air(j));p_air(j+1)=p_air(j); T_air(j+1)=T_air(j); end
% above makes sure rounding errors don't cause v. small expansion
else
p_air(j+1) = r*p_air(j);
end

% calculate temperature of dn moles after compression and compression work
T_air(j+1) = adiabatic_comp(T_air(j), p_air(j+1), p_air(j), k, eta_pol_c);
W_in = dn*adiabatic_work(p_air(j+1), p_air(j), T_air(j), k, eta_pol_c, c_air_molar);
% Calculate the loss of exergy associated with compressing dn moles and sum
% for each stage
Ex_dest_comp = (dn*MM_air)*(heat_cap_air*log(T_air(j+1)/T_air(j))-(R/MM_air)*log(p_air(j+1)/p_air(j)))*T_am;
Comp_Ex_dest(stage)=Comp_Ex_dest(stage)+Ex_dest_comp;

% Record the work asociated with the compression of "dn" at this stage
W_stored(stage) = W_stored(stage) + W_in;    

% through the Packed Bed
T_fluid_in = T_air(j+1);
p_fluid_in = p_air(j+1);
% checking for errors - sometimes snags when temps should be same due to
% rounding error
% if(T_fluid_in<T_am)
% fprintf('air temperature in is too low %f \n', T_fluid_in)
% fprintf('pressure in %f \n', p_air(j+1))
% fprintf('temp pre-comp %f \n', T_air(j))
% fprintf('pressure pre-comp %f \n', p_air(j))
% error('myApp:argChk', 'Stray heat flow')
% end
%Calculate the temperature of dn after it passes through the packed bed

[T_air_out, PB_peb{stage}, PB_fluid{stage}, dn_out, mass_per_dh, p_air(j+2), h_lof, Heat_Ex_loss, ...
    p_bed(stage), moles_in_slice{stage}, Delta_T_fluid{stage}, Delta_n{stage}] = ...
    Reheat_charge(dn, N_dot, T_fluid_in, p_fluid_in, c_air_molar, PB_peb{stage}, PB_fluid{stage}, L_bed(stage), r_bed(stage), dl, c_peb, k_peb, ...
    lambda_ins, ins_th, p_stored, V_stored, T_stored, T_am, ...
    p_bed(stage), moles_in_slice{stage}, Delta_T_fluid{stage}, Delta_n{stage}, dn_out, voidage);

Press_d = p_fluid_in - p_air(j+2);
% % find minimum and maximum values for pressure drop
if(P_drop(stage) == 0); P_drop(stage) = Press_d; end
if(Press_d >= P_drop(stage)); P_drop_max(stage) = Press_d; end
if(Press_d <= P_drop(stage)); P_drop_min(stage) = Press_d; end
% find minumum and maximum values if volumetric heat transfer coefficient
if(h_lof >= h_vol_max); h_vol_max = h_lof; end
if(h_vol_min == 0); h_vol_min = h_lof; end
if(h_lof <= h_vol_min); h_vol_min = h_lof; end

% sum the exergy lost due to thermal power loss from PB
PB_Ex_loss_comp(stage) = PB_Ex_loss_comp(stage) + Heat_Ex_loss;

%The temperature of the fluid after
T_air(j+2)=T_air_out;
% add in a fixed pressure drop term if using one here
% p_air(j+2) = p_air(j+1)-p_drop(stage);

% counting how many moles come out each stage..
n_test(stage) = n_test(stage) + dn_out;

%change dn to represent what comes out
dn = dn_out;

end

% update the amount that actually goes into the store
dn_into_store = dn;

%calculating the new temperature of the air store
Q_stored = n_stored*c_air_molar*T_stored + dn_into_store*c_air_molar*T_air(j+2);
n_stored = n_stored + dn_into_store;

T_stored = Q_stored/(n_stored*c_air_molar);
p_stored = ((n_stored*R*T_stored)/store_volume)/100000;

% recording the pressure and temp at different points
if(p_stored>temp_var)
     moles(temp_var) = n_stored;
     pressure_mat_in(temp_var,:) = p_air;
     temperature_mat_in(temp_var,:) = T_air;
     temp_var = temp_var + 1;
end

    
end

end

%%%%%%%%%%%%%%%% END COMPRESSION PHASE %%%%%%%%%%%%

T_bed_max_check = max(PB_peb{2,1});

% the heat stored in each slice of the PB
Q_contained = zeros(1,NS);
for z = 1:NS
Q_contained(z) = sum((PB_peb{z}-T_am)*mass_per_dh*c_peb);
end
moles(index_max) = n_stored;
moles_max = n_stored;

%include the loss of stratification and temperature loss from the packed
%bed to the ambient

PB_peb_old = PB_peb;
% Accounting for the thermal loss and heat flow through the PB.
for z = 1:NS
[PB_peb{z}, Heat_Ex_loss] = front_collapser_3(PB_peb{z}, k_peb, L_bed(z), r_bed(z), c_peb, dl, ...
    st_time, mass_per_dh, ins_th, lambda_ins, T_am);
PB_Ex_loss_st(z) = PB_Ex_loss_st(z) + Heat_Ex_loss;
end


% logging PB temp profiles
for z = 1:NS
    z1 = num2str(z);    
    name1 = (['Temp_PB_',z1,'_post_comp']);
    eval([name1 '(n_cycle,:) = PB_peb_old{z};']);
    name1 = (['Temp_PB_',z1,'_pre_exp']);
    eval([name1 '(n_cycle,:) = PB_peb{z};']);
end

% if(mod(n_cycle,print_fig)==0)
% %Inspecting the thermal profile of the PB
for z = 1:NS
L=linspace(dl,L_bed(z),Bed_index(z));
figure;
plot(L, PB_peb_old{z}, L, PB_peb{z})
xlabel('length')
ylabel('temperature')
legend('bed temp post comp', 'bed temp pre-exp')
axis([0 max(L) 290 max(PB_peb_old{z})+10])
end
% end

% Flip the order of the bed temperature as during the expansion air flows
% through the cold end first then the hot end
for z = 1:NS
PB_peb{z}=PB_peb{z}(end:-1:1);
PB_fluid{z}=PB_peb{z};
end
clear z;

% maximum temperature and pressure recorded
p_max = p_stored;
T_store_max = T_stored;

% for recording the pressures and temperatures at different stages during
% the expansion
temp_var = ceil(p_max); 
pressure_mat = zeros(temp_var, NS*2+1);
temperature_mat = zeros(temp_var, NS*2+1);
 
%Now calculate temp losses from HP air store
[T_stored] = Store_temp_loss(T_stored, r_store, L_store, ...
    lambda_ins, ins_th, st_time, n_stored, c_air_molar, T_am);

p_stored = (n_stored*R*T_stored/V_stored)/100000;
T_store_end = T_stored; 
for z = 1:NS; Delta_T_fluid{z}(1:Bed_index(z))=0;end

% Set the pressures of the beds equal to the store pressure & other
% expansion pressures
for xx = 1:NS-1; p_bed(xx)=p_test(xx*2); end
p_bed(NS)=p_stored;
% Re-initialise all the PB parameters as otherwise very complicated...
for z = 1:NS
    PB_fluid{z}=PB_peb{z};
    Delta_T_fluid{z}(1:Bed_index(z))=0;
    entry3 = (p_bed(z)*100000*void_vol_slice(z))./(8.31*PB_fluid{z});
    moles_in_slice{z} = entry3;
    clear entry3;
end

%%%%%%%%%%%%%%%% START EXPANSION PHASE %%%%%%%%%%%%

for i = 1:NS
% Making sure store pressure doesn't fall below p_i + P_drop when dn is
% subtracted from the store, 0.005 can be made smaller as dn is made
% smaller and must be increased with dn
if(p_test(2*(NS-i)+1)+0.005+pressure_drop<p_i)
    p_thres = p_i+0.005; % this is for the approach of p_store back to p_i    
else
    p_thres = p_test(2*(NS-i)+1)+0.005+pressure_drop;
end       

while(p_stored>p_thres)

for stage = i:NS    
if(stage>=i)
j = 2*stage - 1;

% if the stage in question is adjacent to the store then the pressure must
% be constantly adjusted as the pressure in the store depletes
if(stage==i)
% Use the the low pressure limit (dn subtracted) as the high pressure for the expansion
p_stored_low = ((n_stored-dn_in)*R*T_stored/V_stored)/100000; 
p_air_ex(j) = p_stored_low; 
% enter the low pressure value
p_air_ex(j+2) = p_test(2*(NS-i)+1);
if(p_air_ex(j)<p_air_ex(j+2)); p_air_ex(j)=p_air_ex(j+2); end % to get a 0 work value
T_air_ex(j)=T_stored;
% the input temp of air into packed beds is store temp
n_stored = n_stored - dn_in;
p_stored = ((n_stored*R*T_stored)/V_stored)/100000;
dn = dn_in;
end
if(stage>i)
    % if not the stage adjacent to store then enter low pressure val
p_air_ex(j+2)=p_test(2*(NS-stage)+1);
end
    
% release dn moles through the subsequent Packed Bed
T_fluid_in = T_air_ex(j);
p_fluid_in = p_air_ex(j);

[T_air_out, PB_peb{NS+1-stage}, PB_fluid{NS+1-stage}, dn_out, mass_per_dh, p_air_ex(j+1), h_lof, Heat_Ex_loss, ...
    p_bed(NS+1-stage), moles_in_slice{NS+1-stage}, Delta_T_fluid{NS+1-stage}] = ...
    Reheat_discharge(dn, N_dot, T_fluid_in, p_fluid_in, c_air_molar, PB_peb{NS+1-stage}, PB_fluid{NS+1-stage}, ...
    L_bed(NS+1-stage), r_bed(NS+1-stage), dl, c_peb, k_peb, lambda_ins, ins_th, p_stored, V_stored, T_stored, T_am, ...
    p_bed(NS+1-stage), moles_in_slice{NS+1-stage}, Delta_T_fluid{NS+1-stage}, voidage);

Press_d = p_fluid_in - p_air(j+2);
% Recording min and max pressure drops
% if(Press_d >= P_drop(NS+stage)); P_drop_max(NS+stage) = Press_d; end
% if(Press_d <= P_drop(NS+stage)); P_drop_min(NS+stage) = Press_d; end
% if(P_drop(stage) == 0); P_drop(stage) = Press_d; end
% % find minumum and maximum values if volumetric heat transfer coefficient
% if(h_lof >= h_vol_max); h_vol_max = h_lof; end
% if(h_vol_min == 0); h_vol_min = h_lof; end
% if(h_lof <= h_vol_min); h_vol_min = h_lof; end

% the moles out is different from the moles in.
dn = dn_out;

% sum the exergy lost due to thermal power loss from PB
PB_Ex_loss_turb(stage) = PB_Ex_loss_turb(stage) + Heat_Ex_loss;

%temperature of the air after
T_air_ex(j+1) = T_air_out;

%option to use a fixed pressure drop
% p_air_ex(j+1) = p_air_ex(j)-p_drop(stage);

% now pass through the turbine
W_out = dn*adiabatic_work_exp(p_air_ex(j+2), p_air_ex(j+1), T_air_ex(j+1), k, eta_pol_t, c_air_molar);

% checking for errors...
test=isnan(W_out);
if(test==1)
% fprintf('%f \n', test)
error('myApp:argChk', 'nan')
end

% calculate new air pressure after expansion
T_air_ex(j+2) = adiabatic_exp(T_air_ex(j+1), p_air_ex(j+2), p_air_ex(j+1), k, eta_pol_t);
W_released(stage) = W_released(stage) + W_out;
Ex_dest_turb = (dn*MM_air)*(heat_cap_air*log(T_air_ex(j+2)/T_air_ex(j+1))-(R/(MM_air))*log(p_air_ex(j+2)/p_air_ex(j+1)))*T_am;
Turb_Ex_dest(stage)=Turb_Ex_dest(stage)+Ex_dest_turb;

end
end

%%% Loss exergy from system exit temperature
Ex_exhaust = dn*MM_air*heat_cap_air*T_am*((T_air_ex(j+2)/T_am)-1-log(T_air_ex(j+2)/T_am));
Exit_loss = Exit_loss + Ex_exhaust;


if(p_stored<temp_var)
     pressure_mat(temp_var,:) = p_air_ex;
     temperature_mat(temp_var,:) = T_air_ex;
     temp_var = temp_var - 1;
end

end
end


%%%%%% END EXPANSION PHASE  %%%%%%%%%%

% %Inspecting the thermal profile of the PB
for z = 1:NS
L=linspace(dl,L_bed(z),Bed_index(z));
figure;
plot(L, PB_peb{z})
xlabel('length')
ylabel('temperature')
legend('bed temp post exp')
axis([0 max(L) 200 500])
end

% Calculating the exergy remaining in the PB
Ex_leftover = cell(A);
for z = 1:NS
    Ex_leftover{z}=zeros(1, Bed_index(z));
end

for z = 1:NS
for i = 1:Bed_index(z)
Ex_leftover{z}(i) = mass_per_dh*c_peb*T_am*((PB_peb{z}(i)/T_am)-1-log(PB_peb{z}(i)/T_am));
end
end
Exergy_left = zeros(1,NS);
for z = 1:NS
    Exergy_left(z) = sum(Ex_leftover{z});
end

Efficiency = sum(W_released)/sum(W_stored);

% cleaning the pressure and temperature matricies so they can be plotted
for i = 1:NS*2+1
for j = 1:ceil(p_max)
if(pressure_mat(j,i) == 0)
    pressure_mat(j,i) = nan;
end
if(temperature_mat(j,i) == 0)
    temperature_mat(j,i) = nan;
end
end
end

% for i = 1:(NS-1)*2
% marker = find(pressure_mat(:,i) == min(pressure_mat(:,i)));
% pressure_mat(marker,i)=nan;
% temperature_mat(marker,i)=nan;
% end
i = 1;
while i<p_i
    pressure_mat_in(i, :) = nan;
    i = i + 1;
end

fprintf('Efficiency estimate = %f \n',-Efficiency*100);
fprintf('Compression work = %f \n',sum(W_stored)/kWh);
fprintf('Work released = %f \n',sum(W_released)/kWh);
fprintf('Exergy destruction compressor = %f \n',sum(Comp_Ex_dest)/kWh);
fprintf('Exergy destruction turbine = %f \n',sum(Turb_Ex_dest)/kWh);
fprintf('Exergy loss packed bed = %f \n',(sum(PB_Ex_loss_comp)+sum(PB_Ex_loss_st)+sum(PB_Ex_loss_turb))/kWh);
fprintf('Exergy left in packed beds = %f \n',sum(Exergy_left)/kWh);
fprintf('Exergy Exit loss = %f \n',Exit_loss/kWh);
fprintf('Exergy Unaccounted = %f \n',(sum(W_stored)+sum(W_released)-sum(PB_Ex_loss_comp)-sum(PB_Ex_loss_st)-sum(PB_Ex_loss_turb)-...
    sum(Comp_Ex_dest)-sum(Turb_Ex_dest)-sum(Exergy_left)-Exit_loss)/kWh);
% Unaccounted exergy is due to pressure losses in the packed beds as well as heat conduction in the packed beds

figure;
x1 = [0 11/2]; y1 = [p_i p_i];
x2 = [0 11/2]; y2 = [p_max p_max];
p2=interp2(log(pressure_mat),meshgrid((1:0.1:5)/2,(1:0.1:80)/2));
surf(p2,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
hold on
p3 = plot(x1, y1, '--k', 'LineWidth', 1.5);
hold on
p4 = plot(x2, y2, '--k', 'LineWidth', 1.5);
xlabel('stage', 'FontSize',12)
ylabel('storage pressure p_s (bar)', 'FontSize',12)
zlabel('ln(pressure)', 'FontSize',12)
title('Illustrating pressure with stage and storage pressure for the 2-stage system during discharge')

% figure;
% rangeY=floor(min(pressure_mat(:,1))):.2:ceil(max(pressure_mat(:,1)));
% rangeX=1:0.1:NS*2+1;
% [X,Y]=meshgrid(rangeX,rangeY);
% Z=interp2(pressure_mat,X,Y,'cubic');

% surf(Z,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
 
% % calculate the max and min charging/discharging power
% % Finding the average power for first second
% n_max= (p_max*100000*store_volume)/(8.31*T_am);
% flowrate = (n_max-n_i)/(charge_time*3600); %moles per sec
% p_1_sec = ((n_i+flowrate)*8.31*T_am)/(100000*store_volume);
% Av_power_1_sec = ((adiabatic_work(p_test(2), p_test(1), T_air(1), k, eta_pol_c, c_air_molar) + ...
%     adiabatic_work(p_test(4), p_test(3), T_air(3), k, eta_pol_c, c_air_molar) + ...
%     adiabatic_work((p_i+p_1_sec)/2, p_test(5), T_air(5), k, eta_pol_c, c_air_molar))*flowrate)/1000;
% % Finding the average power for the last second
% p_last_sec = ((n_max-flowrate)*8.31*T_am)/(100000*store_volume);
% Av_power_last_sec = ((adiabatic_work(p_test(2), p_test(1), T_air(1), k, eta_pol_c, c_air_molar) + ...
%     adiabatic_work(p_test(4), p_test(3), T_air(3), k, eta_pol_c, c_air_molar) + ...
%     adiabatic_work((p_max+p_last_sec)/2, p_test(5), T_air(5), k, eta_pol_c, c_air_molar))*flowrate)/1000;


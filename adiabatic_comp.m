function [T_f] = adiabatic_comp(T_i, p_f, p_i, k, eta_pol)
%Adiabatic temp calculator

%%%% calculate the adjusted polytropic index

T_f = T_i*((p_f/p_i)^((k-1)/(k*eta_pol)));

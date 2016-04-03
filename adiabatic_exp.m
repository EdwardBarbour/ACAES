function [T_f] = adiabatic_exp(T_i, p_f, p_i, k, eta_pol)
%Adiabatic temp calculator

T_f = T_i * (p_f/p_i)^((eta_pol*(k-1))/k);
 
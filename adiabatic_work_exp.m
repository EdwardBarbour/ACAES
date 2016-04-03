function [W] = adiabatic_work_exp(P_f, P_i, T_i, k, eta_pol, c_p)

W = c_p * T_i * ((P_f/P_i)^((eta_pol*(k-1))/k) - 1);

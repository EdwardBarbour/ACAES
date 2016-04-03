function [W] = adiabatic_work(P_f, P_i, T_i, k, eta_pol, c_p)

W = c_p * T_i * ((P_f/P_i)^((k-1)/(eta_pol*k)) - 1);



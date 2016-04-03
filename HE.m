function [T_ho, T_co] = HE(T_hi, E, T_ci)
%Modeling the heat exchanger

dT_max = T_hi - T_ci;

T_ho = T_hi - E * dT_max;

T_co = T_ci + E * dT_max;
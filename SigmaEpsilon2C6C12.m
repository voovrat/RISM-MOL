function [C_o6,C_o12,C_h6,C_h12] = SigmaEpsilon2C6C12(sigma_sl_o, epsilon_sl_o, sigma_sl_h, epsilon_sl_h)

C_o6 = -4*epsilon_sl_o.*sigma_sl_o.^6;
C_o12 = 4*epsilon_sl_o.*sigma_sl_o.^12;

C_h6 = - 4*epsilon_sl_h.*sigma_sl_h.^6;
C_h12 =   4*epsilon_sl_h.*sigma_sl_h.^12;

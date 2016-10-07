function [sigma_new_o, epsilon_new_o, sigma_new_h, epsilon_new_h] = oplsaa_geometric_rules(sigma_in, epsilon_in, sigma_o, epsilon_o, sigma_h, epsilon_h)
% [sigma_new, epsilon_new] = oplsaa_rules(sigma_in, epsilon_in, sigma_o, epsilon_o, sigma_h, epsilon_h);
% oplsaa mixting rule for Solute - O (sigma_o, epsilon_o) and Solute - H (sigma_h, epsilon_h) LJ potential;

sigma_new_o = (sigma_in + sigma_o)/2;
epsilon_new_o = sqrt( epsilon_in.*epsilon_o);

sigma_new_h = ( sigma_in + sigma_h)/2;
epsilon_new_h = sqrt( epsilon_in.*epsilon_h);



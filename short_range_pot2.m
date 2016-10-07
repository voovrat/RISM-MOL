function [EXPONENTA, u_short] = short_range_pot2(LJ_SL_SV, f_s, Bulk_beta, BIG)
% f_s - short -range function for the Ng - procedure.

    u_short =  Bulk_beta * ( LJ_SL_SV + f_s);
    EXPONENTA = exp(-u_short);
    
    I = find( (u_short>BIG)  | isnan(u_short) );
    EXPONENTA(I) = 0;
    u_short(I) = BIG;
    
  
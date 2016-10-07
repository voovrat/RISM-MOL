function [EXPONENTA, u_short] = short_range_pot3(LJ_SL_SV, f_s, Bulk_beta, BIG)
% f_s - short -range function for the Ng - procedure.

%g_treshold = 11125;
%f_split =@(x)( 1./(1+exp(10*(x-log(g_treshold)))) );
%myexp = @(x)( exp(x).*f_split(x) + (x+g_treshold).*(1-f_split(x)) );

    u_short =  Bulk_beta * ( LJ_SL_SV + f_s);
    EXPONENTA = myexp(-u_short);
    
    I = find( (u_short>BIG)  | isnan(u_short) );
    EXPONENTA(I) = 0;
    u_short(I) = BIG;
    
  
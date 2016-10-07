function E=count_err(g0,g1,dR)

%r0 = (1:size(g0,1))*dR;
%r1 = (1:4096)*0.05;
%delim = 0.4/dR;
%f0 = change_grid2(g0,r0,r1,0,0);
%f1 = change_grid2(g1,r0,r1,0,0);
%f0 = g0(1:delim:end,:);
%f1 = g1(1:delim:end,:);


E = mean( sqrt( sum((g0-g1).^2) )*dR );
%E = mean( sqrt( sum((f0-f1).^2) )*0.4 );
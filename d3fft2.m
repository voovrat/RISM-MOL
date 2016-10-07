function y=d3fft2(x,R)
[S,n]=size(x);

dR=R(2)-R(1);
dk=pi/(R(end)+dR);

k = (dk:dk:S*dk)'*ones(1,n);
r = R*ones(1,n);

y=4*pi./k.*dst2(x.*r)*dR;

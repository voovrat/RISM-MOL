function y=d3ifft2(x,R)
[S,n]=size(x);

dR=R(2)-R(1);

dk=pi/(R(end)+dR);

k=(dk:dk:S*dk)'*ones(1,n);
r = R*ones(1,n);

y=1/2/pi^2*dst2(x.*k)./r*dk;
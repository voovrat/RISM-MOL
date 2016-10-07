function Y1=change_grid2(Y0,X0,X1,left_val,right_val,dX)
%
%  written by Volodymyr Sergiievskyi
%  Max Plank Institute for Mathematics in the Sciences
%  2011
%
%
% interpolates function from one grid to another
%
%How to use it:
%
%        Y1=change_grid2(Y0,X0,X1,left_val,right_val,dX)
%
%where Y0 - initial values of the function on grid X0 (column vector/ vectors)
%          X0 - grid of initial function           (column vector/ vectors)
%           X1 - grid of the interpolated function (column vector/ vectors)
%           Y1 - values of the interpolated function  (column vector/ vectors)
%           left_val - if X1 has points, which are smaller than min(X0)-dX, the value of interpolated function in these points will be constant left_val
%          right_val - if X1 has points, which are larger than max(X0)+dX, the value of interpolated function in these points will be constant right_val
%          dX - by default - X0(2)-X0(1), but you may use any value (see left_val, right_val for details) 

[m,n]=size(Y0);

if nargin<6
    dX = X0(2)-X0(1);
end

if nargin==3
    left_val = zeros(n,1);
    right_val = ones(n,1);
end

if length(left_val)==1
    left_val = ones(n,1)*left_val;
end

if length(right_val)==1
    right_val = ones(n,1)*right_val;
end

Y1 = zeros( length(X1) ,n );

for i=1:n

    pp = spline(X0,Y0(:,i));

    x_min0 = min(X0)-dX;
    x_max0 = max(X0)+dX;

    I_INTER = find((X1>=x_min0) & (X1<=x_max0));
    I_LEFT  = find(X1<x_min0);
    I_RIGHT = find(X1>x_max0);

    X_INTER = X1(I_INTER);
    X_LEFT = X1(I_LEFT);
    X_RIGHT = X1(I_RIGHT);

    Y_INTER = ppval(pp,X_INTER);
    Y_LEFT = left_val(i) * ones(size(X_LEFT));
    Y_RIGHT = right_val(i) * ones(size(X_RIGHT));

    Y1( I_LEFT, i ) = Y_LEFT;
    Y1( I_INTER, i ) = Y_INTER;
    Y1( I_RIGHT, i ) = Y_RIGHT;
    
end
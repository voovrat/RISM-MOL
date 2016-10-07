function D = mydist(X)
%
% D=mydist(X). 
%
% lines of X - coordinated of points
%

[m,n]=size(X);

D=zeros(m,m);

for i=1:m
for j=1:m
	d=(X(i,:)-X(j,:));
	D(i,j)=d*d';
end
end

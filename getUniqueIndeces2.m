function [INDEX_MATRIX, VALUES ] = getUniqueIndeces2(D)
%
% D - n x n matix of distances
%
% INDEX_LIST: n*n x 3 array
% each line:  [ i, j, index ] - coordinates of value from D, and it's
% index in values array
%
% value D(i,j) can be found by VALUES(index) = D(i,j) 
%

[m,n]=size(D);

R = reshape(D,m*n,1);

[SR,I]=sort(R);

ROW = (1:m)'*ones(1,n);
COL = ones(m,1)*(1:n);

ROW = reshape(ROW,m*n,1);
COL = reshape(COL,m*n,1);

INDEX_LIST=zeros(m*n,3);

INDEX_LIST(:,1) = ROW(I);
INDEX_LIST(:,2) = COL(I);

INDEX_LIST(1,3) = 1;

val0=SR(1);


VALUES=zeros(m*n,1);
VALUES(1) = val0;

count=1;

for i=2:m*n
    
    if SR(i)~=val0
        val0=SR(i);
        count=count+1;
        
        VALUES(count)=val0;
    end
        
    INDEX_LIST(i,3) = count;
    
end

VALUES=VALUES(1:count);

INDEX_MATRIX = sparse(INDEX_LIST(:,1),INDEX_LIST(:,2),INDEX_LIST(:,3))+zeros(m,n);
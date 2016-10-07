function br = buildRepulsiveBridge_old(C_o12,C_h12,R,W_OH,Bulk_beta)

L=length(C_o12);
N = length(R);

R12 = R.^12;
r12 = R12 * ones(1,L);

C_o12N = ones(N,1)*C_o12';
C_h12N = ones(N,1)*C_h12';
W_OHN = W_OH*ones(1,L);
% repulsive hydrogen and oxygen bridge by Kovalenko w_oh*exp[-Bulk_beta*u_ih]
a=1;

EUST = exp(-Bulk_beta * C_o12N ./ r12 ) - 1;
EUSH = exp(-Bulk_beta * C_h12N ./ r12 ) - 1;

USTB = d3ifft2( d3fft2( EUST, R ).* W_OHN , R);
USTH = d3ifft2( d3fft2( EUSH, R ).* W_OHN , R);

USTB = -USTB./ ( ones(N,1)*USTB(1,:) ) + 1;  
USTH = -USTH./ ( ones(N,1)*USTH(1,:) ) + 1;

br_h = -a * log( abs(USTB ) );
br_h( abs(1-abs(USTB)) < 1e-6) = 0;

I =abs(USTB) < 1e-12;
br_h( I) = Bulk_beta * C_o12N(I) ./ r12(I);

br_o = -a*log( abs(USTH));
br_o( abs(1-abs(USTH)) < 1e-6) = 0;

I = abs(USTH) < 1e-12;
br_o(I) = Bulk_beta * C_h12N(I) ./ r12(I);

br_o( br_o < 1e-12 ) = 0;

br = zeros(N,2*L);

br(:,1:2:end) = br_o / Bulk_beta;
br(:,2:2:end) = br_h / Bulk_beta;

br(isnan(br))=0;

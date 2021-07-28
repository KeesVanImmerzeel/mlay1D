cd "C:/Users/KvI/OneDrive/myR_NOT_on_GitHub/Inspiration/Eendimensionale stroming in meerlagensystemen"

% Analytic steady-state, one-dimensional, leaky multi aquifer model (Matrix functions) 
% Using n connected sections. Left and right most sections run to -inf and +inf resp. 

% input: 
% c, (nLay by nSec) matrix of vertical resistance values of overlaying aquitards 
c = [50, 400, 85/0.075]';
c = [c,c,c,c,c,c,c,c,c,c,c];

% T, (nLay by nSec) matrix of transmissivity values
T = [35*30, 80*30, (55/2)*0.075]';
T = [T,T,T,T,T,T,T,T,T,T,T];

% Q, (nNod by nSec) matrix of nodal injections [L2/T]
Q = zeros(size(c(:,1:end-1)));

% x, (nNod by 1) vector of coordinates of intersection points except +/-inf
x = [-1000,1000,3250,4500,5500,6500,7250,8750,9750,10500];

% X, vector of points where values will be computed 
X = [-2500:10:11000];

% h, (1 by nSec) vector of given heads on top of each sections 
h = [-1.10, -3.85, -1.20, -1.00, -0.80, -0.40, -0.00, 0.40, 0.80, 1.20, 1.60];
% z = [0,-10,-40,-70,-145,-150,-200]';

nLay = length(T(:,1));
nSec = length(T(1,:));

Q = [zeros(nLay, 1), Q, zeros(nLay, 1)];
x = [-inf, x, inf];
Nx = length(x);
H = ones(nLay,1)*h;

xMidSec = 0.5 * (x(1:end-1) + x(2:end));
xMidSec(1) = x(2);
xMidSec(end) = x(end-1);

A = [];
for iSec = 1:nSec 
    a = 1./(T(:,iSec).*c(:,iSec));
    b = 1./(T(:,iSec).*[c(2:nLay,iSec);inf]);
    A(:,:,iSec) = diag(a+b) - diag(a(2:nLay),-1) - diag(b(1:nLay-1),1);
end 

% Write matrix A (3x3x11)
filename = "A.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, A, precision="float64");
fclose (fid);

C = zeros(nLay*(2*(Nx-2)), nLay*(2*(Nx-2)+2));
R = zeros(nLay*(2*(Nx-2)),1);

nNod = nSec-1;
for i = 1:nNod 
    j = i+1; 
    ii = 2*nLay*(i-1)+1; 
    jj = ii+nLay-1; 
    C(ii:jj, ii:jj) = +expm(-(x(j)-xMidSec(i))*sqrtm(A(:,:,i)));
    C(ii:jj, ii+nLay:jj+nLay) = +expm(+(x(j)-xMidSec(i))*sqrtm(A(:,:,i))); 
    C(ii:jj, ii+2*nLay: jj+2*nLay) = -expm(-(x(j)-xMidSec(j))*sqrtm(A(:,:,j))); 
    C(ii:jj, ii+3*nLay: jj+3*nLay) = -expm(+(x(j)-xMidSec(j))*sqrtm(A(:,:,j))); 
    R(ii:jj) = -H(:,i) + H(:,j);
    C(ii+nLay:jj+nLay, ii:jj) = -diag(T(:,i ))*sqrtm(A(:, :,i))*expm(-(x(j)-xMidSec(i))*sqrtm(A(:,:,i)));
    C(ii+nLay:jj+nLay, ii+ nLay:jj+ nLay) = +diag(T(:,i))*sqrtm(A(:,:,i))*expm(+(x(j)-xMidSec(i))*sqrtm(A(:,:,i)));
    C(ii+nLay:jj+nLay, ii+2*nLay:jj+2*nLay) = +diag(T(:,j))*sqrtm(A(:,:,j))*expm(-(x(j)-xMidSec(j))*sqrtm(A(:,:,j)));
    C(ii+nLay:jj+nLay, ii+3*nLay:jj+3*nLay) = -diag(T(:,j))*sqrtm(A(:,:,j))*expm(+(x(j)-xMidSec(j))*sqrtm(A(:,:,j)));
    R(ii+nLay:jj+nLay) = Q(:,j);
end

% Write matrix C (60x66)
filename = "C.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, C, precision="float64");
fclose (fid);

% Write vector R (60x1)
filename = "R.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, R, precision="float64");
fclose (fid);

COEF = C(:,nLay+1:end-nLay)\R;
%Ainv = inv(C(:,nLay+1:end-nLay));
%COEF = [zeros(nLay,1);Ainv * R;zeros(nLay,1)];
COEF = [zeros(nLay,1);COEF;zeros(nLay,1)];

% Write vector COEF (66x1)
filename = "COEF.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, COEF, precision="float64");
fclose (fid);

% output: 
% phi computed heads [H], a (nLay by length(X)) matrix 
% q computed flows [L2/T], a (nLay by length(X)) matrix 
% s [L/T] is the downward positive seepage rate through top of each layer, a nLay by length(X) matrix 
phi = zeros(nLay,length(X));
q = zeros(nLay,length(X));
s = zeros(nLay,length(X)); 

for i=1:length(X)
    iSec = find(X(i) > x(1:end-1) & X(i) <= x(2:end))
    k = 2*nLay*(iSec-1)+1;
    l = k+nLay-1;
    
    sqm = sqrtm(A(:,:,iSec))
    test = expm(-(X(i)-xMidSec(iSec))*sqrtm(A(:,:,iSec)))
    C1 = expm(-(X(i)-xMidSec(iSec))*sqrtm(A(:,:,iSec)))*COEF(k:l); 
    C2 = expm(+(X(i)-xMidSec(iSec))*sqrtm(A(:,:,iSec)))*COEF(k+nLay:l+nLay); 
    C3 = sqrtm(A(:,:,iSec))*C1; 
    C4 = sqrtm(A(:,:,iSec))*C2;  

    phi(:,i) = C1 + C2 + H(:,iSec); 
    q(:,i) = diag(T(:,iSec))*(C3-C4);
    sNet = diag(T(:,iSec))*sqrtm(A(:,:,iSec))*(C3+C4); 

    s(nLay,i) = sNet(nLay);  
    for iLay=nLay-1:-1:1 
        s(iLay,i) = sNet(iLay)+s(iLay+1,i);  
    end 
end 

% Write matrix phi (3x1351)
filename = "phi.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, phi, precision="float64");
fclose (fid);

% Write matrix q (3x1351)
filename = "q.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, q, precision="float64");
fclose (fid);

% Write matrix s (3x1351)
filename = "s.bin";
fid = fopen (filename, "w");
  # Do the actual I/O here…
  i = fwrite (fid, s, precision="float64");
fclose (fid);

plot(X,phi)
plot(X,q(3,:))
plot(X,s(1,:))
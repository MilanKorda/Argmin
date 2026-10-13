% 2D examples of use of the POLYARGMIN function reported in
% D. Henrion, M. Korda, J. B. Lasserre "Polynomial argmin for recovery
% and approximation of multivariate discontinuous functions"

close all
clear all

rng(41421)

% f = @(x1,x2)(x1.^2+x2.^2<=1/4); % disk indicator

f = @(x1,x2)((x1.^2+x2.^2<=1/4)+...
    1/2*((x1-3/4).^2+(x2-3/4).^2<=1/16)+...
    1/2*((x1+3/4).^2+(x2+3/4).^2<=1/16)); % weighted sum of indicators of three disks


% Generate samples
[X1,X2] = meshgrid(linspace(-1,1,25)); % equidistant
Y = double(f(X1,X2));

dx = 6; dy = 6;
P = polyargmin([X1(:)';X2(:)'],Y(:)',dx,dy);
 
% Construct polynomial with GloptiPoly
mpol x 2
mpol y
p = y.^(1:dy) * P * mmon(x,dx);

% Evaluate polynomial on a grid
ymin = min(Y(:)); ymax = max(Y(:));
X1L = linspace(-1,1,200); X2L = linspace(-1,1,200); YL = linspace(ymin,ymax,200);
[X1,X2,YY] = meshgrid(X1L,X2L,YL);
pv = eval(vectorize(p,'X1','X2', 'YY'));

% Compute argmin
[~,I] = min(pv,[],3);
Yargmin = YL(I);

% Display true function and argmin approximation
[X1,X2] = meshgrid(X1L,X2L);
Ytrue = double(f(X1,X2));
mesh(X1,X2,Ytrue); title('True')
figure
mesh(X1,X2,Yargmin); title('Argmin')




% 1D examples of use of the POLYARGMIN function reported in
% D. Henrion, M. Korda, J. B. Lasserre "Polynomial argmin for recovery
% and approximation of multivariate discontinuous functions"

close all
clear all

rng(41421)

% Function to learn

% f = @(xx)(    (2*double(xx > 0.75 | xx < -0.75) - 1) );                        % 2 jumps
% f = @(xx)(    xx.^2.*(2*double(xx > 0.75 | xx < -0.75) - 1) );                 % 2 jumps * x^2
% f = @(xx)(2*double((xx > 0.75 | xx < -0.75) | (xx > -0.25 & xx < 0.25)) - 1);  % 4 jumps
% f = @(xx)(    sin(2*xx).*(2*double(xx > 0.75 | xx < -0.75) - 1) );             % 2jumps * sin
% f = @(x)(1 ./ (1 + 25*x.^2));                                                  % Runge function
% f = @(x)(sqrt(abs(sin(x))));                                                   % sqrt(|sin(x)|)

%
f =  @(x)(     -2/5+2/5*(3./(5-4*cos((x+1)*pi))...                               % 7 discontinuities
               +2*(((x+1)*pi>1)&((x+1)*pi<2.5))...
               -1*(((x+1)*pi>2.5)&((x+1)*pi<3))...
               +1*(((x+1)*pi>3)&((x+1)*pi<4))...
               -1*(((x+1)*pi>4)&((x+1)*pi<4.7))...
               -2*(((x+1)*pi>4.7)&((x+1)*pi<5))...
               -1*(((x+1)*pi>5)&((x+1)*pi<5.4))));
%}

% Generate samples
N = 200; 
X = linspace(-1,1,N);
Y = f(X);

% Construct argmin polynomial
dx = 6; dy = 6;
P = polyargmin(X,Y,dx,dy);

% Construct polynomial with GloptiPoly
mpol x y
p = y.^(1:dy) * P * mmon(x,dx);

% Evaluate polynomial on a grid
ymin = min(Y); ymax = max(Y);
XL = linspace(-1,1,200); YL = linspace(ymin,ymax,200);
[X,Y] = meshgrid(linspace(-1,1,200),YL);
pv = eval(vectorize(p,'X','Y'));

% Compute argmin
[~,I] = min(pv);
Yargmin = YL(I);

% Plot true function and argmin approximation
Ytrue = f(XL);
stairs(XL,Ytrue,'linewidth',5,'color','black'); % True
hold on;
stairs(XL,Yargmin,'linewidth',2,'color','red'); % Argmin
ylim([min([Yargmin]) - 0.1*abs(min([Yargmin])) , max([Yargmin]) + 0.1*abs(max([Yargmin])) ])
legend('True','Argmin')


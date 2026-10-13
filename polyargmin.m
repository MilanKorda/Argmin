% POLYARGMIN   Polynomial argmin approximation of a function
%
% Given a collection of points (column vectors) gathered into
% a matrix X = [x1 x2 .. xN] and the values Y = [F(x1) F(x2) .. F(xN)]
% of a scalar function F at these points, the command
%
%   P = POLYARGMIN(X,Y)
%
% returns a polynomial P such that ARGMIN_Y P(X,Y) approximates F(X).
%
% The polynomial is given as a real matrix P. Each row of P is the 
% coefficient vector of a polynomial Pi(X) in the grevlex ordering, 
% so that Y.^(1:DY) * P * monolist(X, DX) (YALMIP)
% or Y.^(1:DY) * P * mmon(X, DY) (GloptiPoly) returns
% P(X,Y) = P1(X) * Y + P2(X) * Y^2 + .. + PDY(X) * Y^{DY}
%
% P(X,Y) has degree DX in X (default 6) and degree DY in Y (default 6).
% Degrees DX and DY can be specified as additional input arguments.

% Requires: YALMIP and an SDP solver (MOSEK, SeDuMi, SDPT3...)

% Written by D. Henrion, M. Korda, J. B. Lasserre for the experiments
% reported in "Polynomial argmin for recovery and approximation of 
% multivariate discontinuous functions", February 2023.

function P = polyargmin(X,Y,dx,dy)

if nargin < 2
    error('invalid number of input arguments')
end
if nargin < 3
    dx = 6;
end
if nargin < 4
    dy = 6;
end
N = size(X,2);
if size(Y,2) ~= N
    error('Inconsistent column size of first and second input arguments')
end
if size(Y,1) ~= 1
    error('Second input argument must be a row vector')
end

x = sdpvar(size(X,1),1);
y = sdpvar(1,1);

alpha = 1e-2; % curvature parameter

% Construct p parametrized as p = sum_{k=1}^dy hk(x)y^k
p = 0; P = [];
for k = 1:dy
    [h,ch] = polynomial(x,dx);
    p = p + y^k*h;
    P = [P ; ch' ];
end

% Slack polynomial gamma
[gamma,cgamma] = polynomial([x;y],dx);

% SOS multipliers
CS = []; S = [];
for i = 1:N
    [s,cs] = polynomial(y,dy-2);
    CS = [CS ; cs];
    S = [S ; s];
end

% Range of f: Y = {y | g(y) >= 0}
ymin = min(Y); ymax = max(Y);
g = -(y-ymin)*(y-ymax);

% Build constraints and objective
con  = []; obj = 0;
for i = 1:N
    con = [ con, sos(replace( p,x,X(:,i)) - replace(p,[x;y],[X(:,i);Y(i)]) ...
        + replace(gamma,[x;y],[X(:,i);Y(i)]) - g*S(i) - alpha*(y-Y(i))^2 ) ];
    tmp = replace(gamma,[x;y],[X(:,i);Y(i)]);
    con = [con , tmp >= 0];
    obj = obj + tmp;
end
con = [con, sos(S)];

% Solve
solvesos(con,obj,[],[P(:);cgamma;CS])

P = double(P);




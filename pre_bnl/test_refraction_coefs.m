
x=1;
Rf = 1/2.*((x-n.*sqrt(1-((1-x.^2)./n^2))./x+n.*sqrt(1-((1-x.^2)./n^2))).^2+(n.*x-sqrt(1-((1-x.^2)./n^2))./n.*x+sqrt(1-((1-x.^2)./n^2))).^2);

x =   n.*sqrt(1-((1-x.^2)/n^2))

test = @(x) x.*Rf_u;

integral(test,0,1)

% use this to confirm the expression is correct
u = 1;
test = 1/2.*((((u-n.*(sqrt(1-((1-u.^2)./n^2))))./(u+n.*(sqrt(1-((1-u.^2)./n^2))))).^2)+(((n.*u-(sqrt(1-((1-u.^2)./n^2))))./(n.*u+(sqrt(1-((1-u.^2)./n^2))))).^2))

% copy the expression here, multiplied by u, and evaluate the integral
int1 = integral(@(u) u.*1/2.*((((u-n.*(sqrt(1-((1-u.^2)./n^2))))./(u+n.*(sqrt(1-((1-u.^2)./n^2))))).^2)+(((n.*u-(sqrt(1-((1-u.^2)./n^2))))./(n.*u+(sqrt(1-((1-u.^2)./n^2))))).^2)),0,1);
int2 = integral(@(u) u,0,1);

test = int1/int2

% the advantage of integral is that the user does not need to specify du,
% but maybe that is an advantage? 

% repeat using trapz
u = 0:du:1;
int1 = trapz(u,u.*1/2.*((((u-n.*(sqrt(1-((1-u.^2)./n^2))))./(u+n.*(sqrt(1-((1-u.^2)./n^2))))).^2)+(((n.*u-(sqrt(1-((1-u.^2)./n^2))))./(n.*u+(sqrt(1-((1-u.^2)./n^2))))).^2)));
int2 = trapz(u,u);
test = int1/int2;

% try using simplify
syms u
test = simplify(1/2.*((((u-n.*(sqrt(1-((1-u.^2)./n^2))))./(u+n.*(sqrt(1-((1-u.^2)./n^2))))).^2)+(((n.*u-(sqrt(1-((1-u.^2)./n^2))))./(n.*u+(sqrt(1-((1-u.^2)./n^2))))).^2)))

pretty(test)

test(1)
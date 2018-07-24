T = 200;
nvar = 3;
B = 
e = randn(T,nvar);
u = B*e;

% DGP1
% x_{t} - x_{t-1}               = u_{1,t}
% -x_{t} + y_{t} - z_{t}        = u_{2,t} 
% 1/2*x_{t} + 1/2*y_{t} + z_{t} = u_{3,t}
% Two cointegrated vectors: one permanent and two transitory shocks
x = recserar(u(:,1),0,1);
z = -2/3*x-1/3*u(:,2)+2/3*u(:,3);
y = x+z+u(:,2);
DGP1 = [x y z];
clear x y z

% DGP2
% x_{t} - y_{t} - 2*z_{t} = u_{1,t}
% y_{t} - y_{t-1}         = u_{2,t} 
% z_{t} - z_{t-1}         = u_{3,t}
% One cointegrated vector: two permanent and one transitory shocks
y = recserar(u(:,2),0,1);
z = recserar(u(:,3),0,1);
x = y + 2*z + u(:,1);
DGP2 = [x y z];

%DGP King, 

% DGP simple macro model of labor market Hacobson, Vredin and Warne (1997)
% gdp

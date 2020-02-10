function [w_dim] = Wangen2010_variable_EET(h, xc, width, alpha, rho_c, delta_rho, x_dim, g)

load = rho_c*g*h;
Vo = load*width; 
x_dim = transpose(0:1000:max(x_dim)); 
x = abs((x_dim-xc))./alpha;
w = exp(-x).*(cos(x)+sin(x));
w_dim =(((Vo./(2.*delta_rho*g))*w)*-1)./alpha; %dimensionalization




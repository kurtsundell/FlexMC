function [x,w,w_dim_hetenyi] = Hetenyi_variable_EET(k,h_in, xc_in, width_in, alpha_in, rho_c_in, delta_rho_in, x_dim_in,g_in)
%{
for testing
xc_in=xc;
h_in=h_all;
width_in=loadwidth_all;
alpha_in=alpha;
rho_c_in=rho_c;
delta_rho_in=delta_rho;
x_dim_in=x_dim;
g_in=g; 
%}

Vo = rho_c_in*g_in*h_in*width_in; 
x = (0:1000:max(x_dim_in)); 

Ft=transpose(x./alpha_in); %lambda from Hetenyi: x is the dimension of the beam

Fat=xc_in./alpha_in; %lambda-a from Hetenyi: 'a' is the position of the loads
Dzat=(exp(-Fat)).*(cos(Fat)) ; % Functions for origin conditions for loading at 'a' From Hentenyi Ch2 eq h (Hetenyi p. 12)
Czat=exp(-Fat).*[cos(Fat)-sin(Fat)] ; % Functions for origin conditions for loading at 'a' From Hentenyi Ch2 eq h

%preallocate
Wtat = zeros(size(x));%flexure
AbsAzt=zeros(size(x));
a_x=zeros(size(x));
Azt=zeros(size(x));
Bzt=zeros(size(x));

for r=1:size(x,2) ; %for z
    % Functions for origin conditions
    Azt(r)=exp(-Ft(r))*[cos(Ft(r))+sin(Ft(r))] ;%A-lambda-x from Hetenyi Ch 2 equation 23a (Henteny, p. 26)
    Bzt(r)=exp(-Ft(r))*sin(Ft(r)) ; %B-lambda-x from Hetenyi Ch 2 equation 23a (Henteny, p. 26)
    a_x(r)=transpose(abs(x(r)-xc_in)/alpha_in(r));
    AbsAzt(r)=exp(-a_x(r))*[cos(a_x(r))+sin(a_x(r))]; %A-lambda-abs(a-x) from Hetenyi Ch 2 equation 23a (Henteny, p. 26)
    Wtat(r)= [(Czat(r)+(2*Dzat(r)))*Azt(r)-[2*(Czat(r)+Dzat(r))*Bzt(r)]+AbsAzt(r)]; %eq. 23a Hetenyi (1946)
end

done=1;

w=-Wtat.*(Vo./transpose((delta_rho_in)*2*g_in*alpha_in));
w_dim_hetenyi=interp1(x,w,x_dim_in);



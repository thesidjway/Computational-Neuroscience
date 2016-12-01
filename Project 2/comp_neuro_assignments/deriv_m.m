function[m_val]=deriv_m(t,m,V)
phi=1;
alpha_m=-(0.1*phi*(V+35))/(exp(-(V+35)/10)-1);
beta_m=4*phi*exp(-(V+60)/18);

Dm=alpha_m*(1-m)-beta_m*m;
m_val=Dm;


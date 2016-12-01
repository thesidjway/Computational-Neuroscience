function deriv = dydt_hh(t,statevar,Istim)

global F R T RTF ;
global Nao Ko Nai Ki ;
global Cm ;
global GNa GK Gl ENa EK El ;

nvariables = length(statevar) ;
statevarcell = num2cell(statevar) ;

[V,m,h,n] = deal(statevarcell{:}) ;


  %%% compute ionic currents
  gNa = GNa*m^3*h;
  INa = gNa*(V-ENa);

  gK = GK*n^4; 
  IK = gK*(V-EK);

  Il = Gl*(V-El) ;

  Iion = INa + IK + Il ;

  %% compute rate constants to update gates
  alphan = 0.01*(V+50)/(1-exp(-(V+50)/10));
  betan = 0.125*exp(-(V+60)/80);
  
  alpham = 0.1*(V+35)/(1-exp(-(V+35)/10));
  betam = 4*exp(-(V+60)/18) ;
  
  alphah = 0.07*exp(-(V+60)/20);
  betah = 1/(exp(-(V+30)/10)+1);
  
  dm = alpham - (alpham+betam)*m ;
  dh = alphah - (alphah+betah)*h ;
  dn = alphan - (alphan+betan)*n ;
  dV = -(Istim+Iion)/Cm;
  
deriv = ...
  [dV;dm;dh;dn] ;

return



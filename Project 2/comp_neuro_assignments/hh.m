% %% Hodgkin-Huxley model
%    
%    t                   time                    ms
%    V                   membrane potantial      mV
%    INa,IK,Il,Iion      ionic current           uA/cm2
%    Cm                  capacitance             uF/cm2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical constants
global F R T RTF 
F = 96.5;                   % Faraday constant, coulombs/mmol
R = 8.314;                  % gas constant, J/K
T_celsius = 6.3;            % Temperature in celsius
T = 273 + T_celsius ;       % absolute temperature, K 

RTF = R*T/F ;

% default concentrations for squid axon in sea water - mmol/l
global Nao Ko Nai Ki 
Nao = 491 ;
Ko = 20 ;
Nai = 50 ;
Ki = 400 ;

% Cell constant
global Cm 
Cm = 1 ;                            % membrane capacitance, uF/cm^2;

% Maximum channel conductances -- mS/cm^2
global GNa GK Gl ENa EK El 
GNa = 120;
GK = 36;
Gl = 0.3;

% Nernst potentials -- mV
ENa = RTF*log(Nao/Nai);
EK = RTF*log(Ko/Ki);
El =  -49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 2:  Define simulation and stimulus parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tend =  50 ;              % end of simulation, ms

stimdelay = 5 ;
stimdur = 3 ;
stim_amp = 10 ;

stim_start = stimdelay ;
stim_end = stimdelay + stimdur ;

% % % Intervals defined as follows
% % % 1) t=0 zero to beginning of stimulus
% % % 2) beginning to end of stimulus
% % % 3) end of stimulus to end of simulation
simints = 3 ;

intervals(1,:) = [0,stim_start] ;
intervals(2,:) = [stim_start,stim_end] ;
intervals(3,:) = [stim_end,tend] ;

Istim(1) = 0 ;
Istim(2) = stim_amp ;
Istim(3) = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = -60 ;
m = 0 ;
h = 0.6 ;
n = 0.3 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 4:  Loop through and solve model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statevar_i = [V,m,h,n] ;

% % Simulate 60 seconds at rest before stimulus applied
[post,posstatevars] = ode15s(@dydt_hh,[0,60000],statevar_i,[],0) ;
statevar_i = posstatevars(end,:) ;

t = 0 ;
statevars = statevar_i ;
for i=1:simints
  [post,posstatevars] = ode15s(@dydt_hh,intervals(i,:),statevar_i,[],Istim(i)) ;
  t = [t;post(2:end)] ;
  statevars = [statevars;posstatevars(2:end,:)] ;
  statevar_i = posstatevars(end,:) ;
end

outputcell = num2cell(statevars,1) ;

[V,m,h,n] = deal(outputcell{:}) ;

gNa = GNa*m.^3.*h;
INa = gNa.*(V-ENa);

gK = GK*n.^4; 
IK = gK.*(V-EK);

Il = Gl*(V-El) ;
Iion = INa + IK + Il ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Step 5:  Plot or write output to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1)
plot(t,V)
set(gca,'TickDir','Out')
xlabel('time (ms)')
ylabel('V_m (mV)')

subplot(2,2,2)
hold on
plot(t,INa,'b')
plot(t,IK,'r')
plot(t,Il,'g')
set(gca,'TickDir','Out')
legend('I_N_a','I_K','I_L')
xlabel('time (ms)')
ylabel('Ionic current')
hold off

subplot(2,2,3)
hold on
plot(t,gNa,'b')
plot(t,gK,'r')
set(gca,'TickDir','Out')
legend('g_N_a','g_K')
xlabel('time (ms)')
ylabel('Conductances')

subplot(2,2,4)
hold on
plot(t,m,'b')
plot(t,h,'r')
plot(t,n,'g')
set(gca,'TickDir','Out')
legend('m','h','n')
xlabel('time (ms)')
ylabel('Gating variables')



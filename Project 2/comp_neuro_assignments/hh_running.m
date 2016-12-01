global E_Na E_K E_l T_celcius T Cm Q

E_Na=55; %mV
E_K=-72; %mV
E_l=-49; %mV
T_celcius=6.3; 
T=273+6.3;%kelvin
Cm=1; %microFarad/cm^2;
Q=3; 

global G_Na G_K G_l 
G_Na=120; %milliSimens/cm^2;
G_K=36; %milliSimens/cm^2;
G_l=0.3; %milliSimens/cm^2;
I_ext=[0 0 0]; %uA/cm^2

t_tot=50;
curr_start=5;
curr_dur=40; % increase this to generate multiple action potentials
curr_stop=curr_start+curr_dur;

V=-60; %mV
m=0;
h=0.6;
n=0.3;
dynamic_variables=[V m h n];


init_states=[V m h n];
time=0;

[t,vars]=ode15s(@(t,vars) deriv_hh(t,vars,I_ext(1)),[0 curr_start],init_states);
dynamic_variables=[dynamic_variables; vars];
time=[time;t];
[t,vars]=ode15s(@(t,vars) deriv_hh(t,vars,I_ext(2)),[curr_start curr_stop],dynamic_variables(end,:));
dynamic_variables=[dynamic_variables; vars];
time=[time;t];
[t,vars]=ode15s(@(t,vars) deriv_hh(t,vars,I_ext(3)),[curr_stop 50],dynamic_variables(end,:));
dynamic_variables=[dynamic_variables; vars];
time=[time;t];

figure
subplot(2,1,1)
plot(time,dynamic_variables(:,1));
hold on
plot((curr_start:0.01:curr_stop),I_ext(2),'r','LineWidth',5);
grid on
xlabel('Time ms')
ylabel('membrane voltage mV')
title('action potential at 10uA')

subplot(2,1,2)
hold on
plot(time,dynamic_variables(:,2),'b');
plot(time,dynamic_variables(:,3),'r');
plot(time,dynamic_variables(:,4),'g');
set(gca,'TickDir','Out');
legend('m','h','n');
grid on
xlabel('time ms')
ylabel('probabilities')











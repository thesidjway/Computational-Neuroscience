mu=[0.1 1 100];
tic
for i=1:length(mu)
    [t,vars]=ode45(@(t,vars) deriv_vdp(t,vars,mu(i)),[0 500],[1,0]);
    figure
    plot(vars(:,1),vars(:,2));
    title(mu(i));
end
toc

tic
for i=1:length(mu)
    [t,vars]=ode15s(@(t,vars) deriv_vdp(t,vars,mu(i)),[0 500],[1,0]);
    figure
    plot(vars(:,1),vars(:,2));
    title(mu(i));
end
toc


%Question no 1: Gaussian Estimation
k=zeros(100);
j=zeros(100);

for i = -49:50
    if(i>0)
    for l=1:20000-i
        j(i+49)=j(i+49)+(Stimulus(l)*Stimulus(l+i));
    end
    j(i+49)=j(i+49)/(20000-i);
    end
    if (i<0)
    for l=-1*i+1:20000
        j(i+50)=j(i+50)+(Stimulus(l)*Stimulus(l+i));
    end
    j(i+50)=j(i+50)/(20000+1*i-1);
    end
end

t=linspace(-50,50,100);
figure(1)
plot(t,j);
psth=zeros(4,20000);
%Question no 2: PSTH Evaluation
for i = 1:4
    for j = 1:50
        [m,n]=size(All_Spike_Times{i,j});
        for p = 1:n
            temp=ceil(All_Spike_Times{i,j}(p)*1000);
            psth(i,temp)=psth(i,temp)+1;
        end
    end
end
figure(2)
%PSTH for 4 neurons
    ax1=subplot(4,1,1)
    plot(1000*psth(1,1:20000),'r');
    xlabel(ax1,'time (ms)');
    ylabel(ax1,'r(t)');
    ax2=subplot(4,1,2)
    plot(1000*psth(2,1:20000));
    xlabel(ax2,'time (ms)');
    ylabel(ax2,'r(t)');
    ax3=subplot(4,1,3)
    plot(1000*psth(3,1:20000),'r');
    xlabel(ax3,'time (ms)');
    ylabel(ax3,'r(t)');
    ax4=subplot(4,1,4)
    plot(1000*psth(4,1:20000));
    xlabel(ax4,'time (ms)');
    ylabel(ax4,'r(t)');







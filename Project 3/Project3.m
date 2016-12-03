%%================================================================%%
%%  Computational Neuroscience (EC60007) Project 3
%%
%%  Siddharth S Jha,  14EE30022
%%  Adarsh Mukesh,    13BT3
%%  
%%================================================================%%
%%  Question no 1: Gaussian Estimation
k=zeros(100);
j=zeros(100);

for i = -49:50
    if(i>0)
    for l=1:20000-i
        j(i+50)=j(i+50)+(Stimulus(l)*Stimulus(l+i));
    end
    j(i+50)=j(i+50)/(20000-i);
    end
    if (i<=0)
    for l=-1*i+1:20000
        j(i+50)=j(i+50)+(Stimulus(l)*Stimulus(l+i));
    end
    j(i+50)=j(i+50)/(20000+1*i-1);
    end
end

t=linspace(-50,50,100);
figure(1)
plot(t,j);

%%================================================================%%
%%  Question no 2: PSTH Evaluation

psth=zeros(4,20000);
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
    ax1=subplot(4,1,1);
    plot(1000*psth(1,1:20000),'r');
    xlabel(ax1,'time (ms)');
    ylabel(ax1,'r(t)');
    ax2=subplot(4,1,2);
    plot(1000*psth(2,1:20000));
    xlabel(ax2,'time (ms)');
    ylabel(ax2,'r(t)');
    ax3=subplot(4,1,3);
    plot(1000*psth(3,1:20000),'r');
    xlabel(ax3,'time (ms)');
    ylabel(ax3,'r(t)');
    ax4=subplot(4,1,4);
    plot(1000*psth(4,1:20000));
    xlabel(ax4,'time (ms)');
    ylabel(ax4,'r(t)');
    
%%================================================================%%
%%  Question no 3: Poisson or Non-Poisson

bintimes=[10, 20, 50, 100, 200, 500];

for i = 1:6
    noofbin(i)=20000/bintimes(i);
    binfreq(i)=1000/bintimes(i);
end

spikes=zeros(4,6,2000);
for neur_no=1:4
    for freq_no=1:6
        for iter = 1:20000
            poissonmat{neur_no,freq_no,iter}=0;
        end
    end
end

varmat=[];
meanmat=[];

for neur_no=1:4
    for freq_no=1:6
        for j = 1:50
        [m,n]=size(All_Spike_Times{neur_no,j});
            for p = 1:n
                temp=ceil(All_Spike_Times{neur_no,j}(p)*binfreq(freq_no));
                spikes(neur_no,freq_no,temp)=spikes(neur_no,freq_no,temp)+1;
            end
%             poissonmat{neur_no,freq_no}=spikes(neur_no,freq_no,1:noofbin(freq_no));
            for o=1:noofbin(freq_no)
                for b=bintimes(freq_no)*o-bintimes(freq_no)+1:bintimes(freq_no)*o
                    poissonmat{neur_no,freq_no,b}=spikes(neur_no,freq_no,o);
                end
%                     poissonmat{neur_no,freq_no,10*o}=spikes(neur_no,freq_no,o);
            end
            spikes=zeros(4,6,20000);
            varmat(neur_no,freq_no,j)=  var([poissonmat{neur_no,freq_no,:}]);
            meanmat(neur_no,freq_no,j)= mean([poissonmat{neur_no,freq_no,:}]);
        end
    end
end


figure(3) %for 10ms
subplot(2,2,1)
scatter(varmat(1,1,:),meanmat(1,1,:))
subplot(2,2,2)
scatter(varmat(2,1,:),meanmat(2,1,:))
subplot(2,2,3)
scatter(varmat(3,1,:),meanmat(3,1,:))
subplot(2,2,4)
scatter(varmat(4,1,:),meanmat(4,1,:))

figure(4) %for 20ms
subplot(2,2,1)
scatter(varmat(1,2,:),meanmat(1,2,:))
subplot(2,2,2)
scatter(varmat(2,2,:),meanmat(2,2,:))
subplot(2,2,3)
scatter(varmat(3,2,:),meanmat(3,2,:))
subplot(2,2,4)
scatter(varmat(4,2,:),meanmat(4,2,:))

figure(5) %for 50ms
subplot(2,2,1)
scatter(varmat(1,3,:),meanmat(1,3,:))
subplot(2,2,2)
scatter(varmat(2,3,:),meanmat(2,3,:))
subplot(2,2,3)
scatter(varmat(3,3,:),meanmat(3,3,:))
subplot(2,2,4)
scatter(varmat(4,3,:),meanmat(4,3,:))

figure(6) %for 100ms
subplot(2,2,1)
scatter(varmat(1,4,:),meanmat(1,4,:))
subplot(2,2,2)
scatter(varmat(2,4,:),meanmat(2,4,:))
subplot(2,2,3)
scatter(varmat(3,4,:),meanmat(3,4,:))
subplot(2,2,4)
scatter(varmat(4,4,:),meanmat(4,4,:))

figure(7) %for 200ms
subplot(2,2,1)
scatter(varmat(1,5,:),meanmat(1,5,:))
subplot(2,2,2)
scatter(varmat(2,5,:),meanmat(2,5,:))
subplot(2,2,3)
scatter(varmat(3,5,:),meanmat(3,5,:))
subplot(2,2,4)
scatter(varmat(4,5,:),meanmat(4,5,:))


figure(8) %for 500ms
subplot(2,2,1)
scatter(varmat(1,6,:),meanmat(1,6,:))
subplot(2,2,2)
scatter(varmat(2,6,:),meanmat(2,6,:))
subplot(2,2,3)
scatter(varmat(3,6,:),meanmat(3,6,:))
subplot(2,2,4)
scatter(varmat(4,6,:),meanmat(4,6,:))


%%================================================================%%
%%  Question no 4: Spike Triggered Average
averagemat=zeros(4,100);
staspikes=zeros(4,1,200);
Stim15=Stimulus(1:15000);

mean_rate=zeros(4);
for neur_no=1:4
    n_spikes=0;
    for j=1:50
        v=[];
        v=All_Spike_Times{neur_no,j}<15;
        newspike15=All_Spike_Times{neur_no,j}(v);
        [m,n]=size(newspike15);
        n_spikes=n_spikes+n;
        for spk_no=1:n
            for tim=1:100
                averagemat(neur_no,tim)= averagemat(neur_no,tim)+Stim15(max(1,ceil(1000*(newspike15(spk_no)))-tim));
            end
        end
    end
    mean_rate(neur_no)=n_spikes/750;
    for tim=1:100
        averagemat(neur_no,tim)= averagemat(neur_no,tim)/n_spikes;
    end
end


% plot(averagemat(1,:));
% ylim([-0.2 0.2]);
% subplot(2,2,2)
% plot(averagemat(2,:));
% ylim([-0.2 0.2]);
% subplot(2,2,3)
% plot(averagemat(3,:));
% ylim([-0.2 0.2]);
% subplot(2,2,4)
% plot(averagemat(4,:));
% ylim([-0.2 0.2]);

corr_new=zeros(101);
for i = 0:100
    for l=1:20000-i
        corr_new(i+1)=corr_new(i+1)+(Stimulus(l)*Stimulus(l+i));
    end
    corr_new(i+1)=corr_new(i+1)/(20000-i);
end
figure(4)
plot(corr_new);
corr_mat=zeros(100,100);

for row=1:100
    for col=1:100
        corr_mat(row,col)=corr_new(abs(row-col)+1);
    end
end
cssinv=inv(corr_mat);
ratenew=mean_rate(:,1);
corrected_averagemat=cssinv*transpose(averagemat);

for i =1:4
    corrected_averagemat(:,i)=corrected_averagemat(:,i)*ratenew(i)
end


figure(5)
subplot(2,2,1)
plot(corrected_averagemat(:,1));
subplot(2,2,2)
plot(corrected_averagemat(:,2));
subplot(2,2,3)
plot(corrected_averagemat(:,3));
subplot(2,2,4)
plot(corrected_averagemat(:,4));

%%================================================================%%
% %%  Question no 5: Output nonlinearity
% 
% diagvar=zeros(100,100);
% for i=1:100
%     for j=1:100
%         if(i==j)
%             diagvar(i,i)=0.33;
%         end
%     end
% end
% 
% cssinv=inv(diagvar);
% 
% corrected_averagemat=cssinv*transpose(averagemat);%*mean_rate(:,1)
% 
% for i =1:4
%     corrected_averagemat(:,i)=corrected_averagemat(:,i)*ratenew(i)
% end
% figure(1)
% subplot(2,2,1)
% plot(corrected_averagemat(:,1));
% subplot(2,2,2)
% plot(corrected_averagemat(:,2));
% subplot(2,2,3)
% plot(corrected_averagemat(:,3));
% subplot(2,2,4)
% plot(corrected_averagemat(:,4));
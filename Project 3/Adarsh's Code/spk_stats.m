%% load and autocorr
load('CN_Project3_2016.mat');
trn_stm=Stimulus(1,1:15000);
test_stm=Stimulus(1,15001:end);
for i=1:4
    for j=1:50
        trn_spk{i,j}=[];
        test_spk{i,j}=[];
        a=All_Spike_Times{i,j};
        for k=1:length(a)
            if a(1,k)<15
                trn_spk{i,j}=[trn_spk{i,j},a(1,k)];
            else
                test_spk{i,j}=[test_spk{i,j},a(1,k)-15];
            end
        end
    end
end
figure
autocorr(Stimulus,50);
%% PSTH
figure
binsz=0.05;
trn_psth=zeros(4,15/binsz);
for n_n=1:4
    for i=1:50
        a=trn_spk{n_n,i};
        for j=1:length(a)
            trn_psth(n_n,fix(a(j)/binsz)+1)=trn_psth(n_n,fix(a(j)/binsz)+1)+1;
        end
    end
    trn_psth(n_n,:)=trn_psth(n_n,:)./(50*binsz);
    subplot(2,2,n_n);
    bar(trn_psth(n_n,:));
    str=strcat('neuron ',num2str(n_n),' train psth');
    title(str);
end

binsz=0.001;
test_psth=zeros(4,15/0.001);
figure
for n_n=1:4
    for i=1:50
        a=test_spk{n_n,i};
        for j=1:length(a)
            test_psth(n_n,fix(a(j)/binsz)+1)=test_psth(n_n,fix(a(j)/binsz)+1)+1;
        end
    end
    test_psth(n_n,:)=test_psth(n_n,:)./0.05;
    subplot(2,2,n_n);
    bar(test_psth(n_n,:));
    str=strcat('neuron ',num2str(n_n),' test psth');
    title(str);
end

%% Poisson or non poisson

binsz=[0.01,0.02,0.05,0.1,0.2,0.5];
for i=1:length(binsz)
    figure
    for n_n=1:4
        poss{n_n,i}=zeros(50,20/binsz(i));
        for j=1:50
            a=All_Spike_Times{n_n,j};
            for ii=1:length(a)
                poss{n_n,i}(j,fix(a(ii)/binsz(i))+1)=poss{n_n,i}(j,fix(a(ii)/binsz(i))+1)+1;
            end
        end
        poss{n_n,i}=poss{n_n,i}./binsz(i);
        subplot(2,2,n_n)
        scatter([1:1:50],mean(poss{n_n,i},2),5);
        hold on
        scatter([1:1:50],var(poss{n_n,i},0,2),5);
        str=strcat('neuron ',num2str(n_n));
        title(str);
    end
end

%% STA
stim=[zeros(1,99), trn_stm];
sta=zeros(4,100);
figure
for n_n=1:4
    spk=zeros(50,100);
    filt=zeros(1,100);
    k=0;
    for i=1:50
        a=trn_spk{n_n,i};
        for j=1:length(a)
            if fix(a(j)/0.001)>0
                filt=filt+stim(fix(a(j)/0.001):fix(a(j)/0.001)+99);
                k=k+1;
            end
        end
        
    end
    
    sta(n_n,:)=filt./k;
    subplot(2,2,n_n)
    plot(sta(n_n,:))
    str=strcat('STA neuron',num2str(n_n));
    title(str);
end
x=[];
for i=1:15000
    x=[x;stim(i:i+99)];
end
css=x'*x;
css=css./15000;
figure
for n_n=1:4
    sta_w(n_n,:)=(css^-1)*sta(n_n,:)';
    subplot(2,2,n_n)
    plot(sta_w(n_n,:));
    str=strcat('whitened STA neuron',num2str(n_n));
    title(str);
end

%% non-linearity
figure
for n_n=1:4
    y_t(n_n,:)=conv(sta(n_n,:),trn_stm);
    subplot(2,2,n_n);
    scatter(mean(reshape(y_t(n_n,1:15000),50,[]),1),trn_psth(n_n,:)*10,1);
    str=strcat('neuron',num2str(n_n),' nonlinearity');
    title(str);
end


            
        
        
        




    
   
    


            
    
    
 



    
            
            




        
        
        

           
        
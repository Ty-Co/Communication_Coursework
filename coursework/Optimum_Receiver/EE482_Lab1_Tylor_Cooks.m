%% PART 1: DECISION BOUNDRY AND DECISION RULE
% Define parameters

n = 500; % number of samples
s_ij = [1 .5; -1 1]; % signal projected
E1 = 1.25; %Energy signal 1 
E2 = 2;    %Energy power signal 2
data = randi([0,1],[1,n]);
%%
SNRs = [-3 0 3 6 15];
i = 1;
for snr_db = SNRs
    snr = 10^(snr_db/10);
    sigma1 = (E1/snr);
    sigma2 = (E2/snr);
    
    w1 = normrnd(0,sigma1,1,n);
    w2 = normrnd(0,sigma2,1,n);
 
    r11 = (s_ij(1,1)+ w1);
    r12 = (s_ij(1,2)+ w2);
    r21 = (s_ij(2,1)+ w1); 
    r22 = (s_ij(2,2)+ w2);
    
    %if the intended signal is 1, use this
    d1 = (r21 - s_ij(1,1)).^2 + (r22 - s_ij(1,2)).^2 >=...
         (r21 - s_ij(2,1)).^2 + (r22 - s_ij(2,2)).^2;
     
    %if the intended signal is 0, use this
    d0 = (r11 - s_ij(1,1)).^2 + (r12 - s_ij(1,2)).^2 >=...
         (r11 - s_ij(2,1)).^2 + (r12 - s_ij(2,2)).^2;
     
    r = zeros(2,n);
    figure(1)
    subplot(1,5,i) 
    for v = 1:n
      if data(v) == 1
        r(1,v) = (r21(v));
        r(2,v) = (r22(v));
        if d1(v) == data(v)
            scatter(r21(v),r22(v),'g')
            
        else
            scatter(r21(v),r22(v),'r')
        end
      else  
        r(1,v) = (r11(v));
        r(2,v) = (r12(v));
        if d0(v) == data(v)
            scatter(r11(v),r12(v),'b')
        else 
            scatter(r11(v),r12(v),'o')
        end
      end
      hold on
    end
    title(['r: ', num2str(SNRs(i)), ' dB']);
    hold off
    

    if snr_db == 0
        figure(2)
        subplot(1,2,1)
        histogram(r(1,:),100)
        title(['pdf \phi1: ', num2str(SNRs(i)), ' dB']);
        
        subplot(1,2,2)
        histogram(r(2,:),100)
        title(['pdf \phi2: ', num2str(SNRs(i)), ' dB']);
    end

    if snr_db == 15
        figure(3)
        subplot(1,2,1)
        histogram(r(1,:),100)
        title(['pdf \phi1: ', num2str(SNRs(i)), ' dB']);
        
        subplot(1,2,2)
        histogram(r(2,:),100)
        title(['pdf \phi2: ', num2str(SNRs(i)), ' dB']);
    end
    
    
    
    i = i+1;
end
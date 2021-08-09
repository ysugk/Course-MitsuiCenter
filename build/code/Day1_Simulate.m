%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_Simulate.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement baseline simulation')
disp(' ')
tic

%%%%%%%%%% simulate the exogenous process m
disp('Simulating exogenous process')


% draw the U(0,1) random shocks for simulation of exogenous m
firmushocks = rand(firmnum,Ttot);

%form cumulative probability bins for comparison
pr_thresh_z = cumsum(pr_mat_z,2);

%for each firm, do the exogenous simulation of z
zsimpos = zeros(firmnum,Ttot);

%initialize each firm's m to specified value
t=1;
zsimpos(:,t) = zinit;

%loop over firms and periods
for firmct=1:firmnum
    for t=2:Ttot
  
      %what was the last value of z?
      zctmin1 = zsimpos(firmct,t-1);
      
      %what is today's uniform shock?
      tshockval = firmushocks(firmct,t);
      
      %compare the shock to the appropriate thresholds
      lessvec = find(tshockval<=pr_thresh_z(zctmin1,:));
      zct = lessvec(1);
      
      zsimpos(firmct,t) = zct;
      
    end
end

zsim = zgrid(zsimpos);

%%%%%%%%%% simulate the endogenous pricing process
disp('Simulating endogenous process')
ksimpos = zeros(firmnum,Ttot);

t=1;
ksimpos(:,t) = kinit;

for firmct=1:firmnum
   for t=2:Ttot
    
      %extract states
      kmin1ct=ksimpos(firmct,t-1);
      zct = zsimpos(firmct,t-1);
      
      %store policy
      ksimpos(firmct,t) = polind(kmin1ct,zct);
      
   end
end

ksim = kgrid(ksimpos);

%%%%%%%%%% throw away burn-in period
ergsamp = (Terg+1):(Ttot-1);
zsimpos = zsimpos(:,ergsamp);
zsim = zsim(:,ergsamp);
ksimpos = ksimpos(:,ergsamp);
ksim = ksim(:,ergsamp);

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick simulation check: k lb, min(ksim), max(ksim), k ub')
disp(num2str([kgrid(1) min(ksim(:)) max(ksim(:)) kgrid(knum)]))
disp(' ')


%%%%%%%%%% plot some sample paths and distributions
figure;
plot((1:Tsim),ksim(1,:),'b',...
    (1:Tsim),ksim(2,:),'g',...
    (1:Tsim),ksim(3,:),'r',...
    'LineWidth',lwidnum)
xlabel('Time t')
title('Firm Capital')
ylabel('k_t')
axis([1 Tsim min(ksim(:))*0.95 max(ksim(:))*1.05])
set(gca,'FontSize',fsizenum)

figure;
plot((1:Tsim),zsim(1,:),'b',...
    (1:Tsim),zsim(2,:),'g',...
    (1:Tsim),zsim(3,:),'r',...
    'LineWidth',lwidnum)
xlabel('Time t')
title('Firm Profitability')
ylabel('z_t')
axis([1 Tsim min(zsim(:))*0.95 max(zsim(:))*1.05])
set(gca,'FontSize',fsizenum)

figure;
hist(zsim(:))
title('Simulated Distribution: Firm Profitability')
ylabel('Simulated Frequency')
xlabel('m')
set(gca,'FontSize',fsizenum)

figure;
hist(ksim(:))
title('Simulated Distribution: Firm Capital')
ylabel('Simulated Frequency')
xlabel('k')
set(gca,'FontSize',fsizenum)

figure;
scatter(zsim(:),ksim(:),'b')
title('Simulated Scatter Plot')
ylabel('Firm Capital k')
xlabel('Firm Profitability z')
set(gca,'FontSize',fsizenum)

klag = ksim(:,1:(end-1));
know = ksim(:,2:end);
figure;
scatter(klag(:),know(:),'b'); hold on;
plot(kgrid,kgrid,'k','LineWidth',lwidnum); 
axis([min(ksim(:))*0.95 max(ksim(:))*1.05 min(ksim(:))*0.95 max(ksim(:))*1.05])
title('Simulated Scatter Plot')
ylabel('Firm Capital k')
xlabel('Firm Lagged Capital k_{-1}')
set(gca,'FontSize',fsizenum)

toc
disp(' ')

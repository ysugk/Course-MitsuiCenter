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

toc
disp(' ')
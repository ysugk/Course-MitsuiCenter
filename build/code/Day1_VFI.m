%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_VFI.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement baseline VFI solution')
disp(' ')
tic

%%%%%%%%%% do the VFI loop

%initialize guess for value function
Vold = zeros(statenum,1);
polindold = zeros(statenum,1);
polold = polindold;

for vfct=1:maxit
    
    %based on guess for value function, form continuation value EVmat 
    %with (i,j) entry equal to the EV(') at state i given policy j
    EVmat = kron((reshape(Vold,knum,znum)*(pr_mat_z'))',ones(knum,1));
   
    %form RHS matrix with static and continutation payoff for state i and 
    %policy j
    RHSmat = Emat + beta*EVmat;
    
    %actually do the optimization
    [V,polind] = max(RHSmat,[],2);
    pol = kgrid(polind);
    
    %now, compute the error in the value function from one iteration to the next
    solerr = max(abs(V(:)-Vold(:)));
    polerr = max(abs(pol(:)-polold(:)));
    
    %display some diagnostics
    if (mod(vfct,5)==1) 
        disp(['For iter ' num2str(vfct) ', VF error = ' num2str(solerr)])
    end
    
    %if VF has converged, break loop
    if (solerr<soltol)
        break
    end
    
    %otherwise, update VF and continue
    Vold = V;
    polindold=polind;
    polold = pol;
    
end

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick policy check: k lb, min(pol), max(pol), k ub')
disp(num2str([kgrid(1) min(pol) max(pol) kgrid(knum)]))
disp(' ')

%%%%%%%%%% some simple plots
V = reshape(V,knum,znum);
pol = reshape(pol,knum,znum);
polind = reshape(polind,knum,znum);

figure;
plot(kgrid,V(:,1),'b',...
    kgrid,V(:,floor(znum/2)),'g',...
    kgrid,V(:,znum),'r',...
    'LineWidth',lwidnum)
xlabel('Old Capital k_{-1}')
ylabel('Firm Value V(k_{-1},m)')
axis([kgrid(1) kgrid(knum) min(V(:))*0.95 max(V(:))*1.05])
set(gca,'FontSize',fsizenum)
legend('Low z','Medium z','High z','FontSize',fsizenum,'Location','northwest')
legend boxoff

p = gcf;
exportgraphics(p, 'build/output/figs/Day1-fig1.pdf');
clf;

figure
plot(kgrid,pol(:,2),'b',...
    kgrid,pol(:,floor(znum/2)),'g',...
    kgrid,pol(:,znum),'r',...
    kgrid,kgrid,'k',...
    'LineWidth',lwidnum)
xlabel('Old Capital k_{-1}')
ylabel('New Capital k(k_{-1},z)')
axis([kgrid(1) kgrid(knum) kgrid(1) kgrid(knum)])
set(gca,'FontSize',fsizenum)
legend('Low z','Medium z','High z','Identity','FontSize',fsizenum,'Location','northwest')
legend boxoff

p = gcf;
exportgraphics(p, 'build/output/figs/Day1-fig2.pdf');
clf;

toc
disp(' ')

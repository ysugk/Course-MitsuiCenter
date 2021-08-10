%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Project
% Day1_VFI_MPQ.m
% Yongseok Kim - Indiana University
% 2021 Summer Summer School on Structural Estimation in Corporate Finance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement VFI solution with MQP bounds')
disp(' ')
tic

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
    
    %otherwise, compute MQP bounds ot use in update of VF and continue
    
    minshift = (beta/(1-beta))*min(V(:)-Vold(:));
    maxshift = (beta/(1-beta))*max(V(:)-Vold(:));
    Vold = V + (maxshift+minshift)/2;
    polindold=polind;
    polold = pol;
    
end

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick policy check: k lb, min(pol), max(pol), k ub')
disp(num2str([kgrid(1) min(pol) max(pol) kgrid(knum)]))
disp(' ')

V = reshape(V,knum,znum);
pol = reshape(pol,knum,znum);
polind = reshape(polind,knum,znum);

toc
disp(' ')
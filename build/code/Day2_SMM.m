function Q = Day2_SMM(params)
    alpha = params(1);
    delta = params(2);
    
    disp('Current parameter estimates')
    disp(num2str(params))

    global beta rho sig lambda knum znum znstdev statenum ...
    soltol maxit ...
    firmnum Terg Tsim Ttot zinit kinit ...
    datamom nummom

    run('Day2_solve_model');
    run('Day2_simulate');
    
    profsim = zsim .* (ksim .^alpha) ./ ksim;
    klag = ksim(:,1:(end-1));
    know = ksim(:, 2:end);
    invsim = (know - (1-delta)*klag) ./ klag;
    
    simmom = Day2_compute_moments(profsim, invsim);
    
    W = Day2_weight_matrix(nummom);
    
    momerr = datamom - simmom;
    
    Q = momerr * W * momerr';

end
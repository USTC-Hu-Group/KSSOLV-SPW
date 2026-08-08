function opt_report(opt)
fprintf('----------------------options----------------------\n')
fprintf('eig_method: %s.\n',opt.eigmethod);

fprintf('scf_tol Opt:\n')
fprintf('scftol_fine: %2.2e.\n',opt.scftol);
fprintf('scftol_coarse: %2.2e.\n',opt.scftol_coarse);

fprintf('cgtol Opt:\n')

fprintf('cgtol_fine: %2.2e.\n',opt.cgtol);
fprintf('cgtol_coarse:  %2.2e.\n',opt.cgtol_c);

fprintf('maxcgiter Opt:\n')

fprintf('cgiter_fine: %02d.\n',opt.maxcgiter);
fprintf('cgiter_coarse:  %02d.\n',opt.maxcgiter_c);

end
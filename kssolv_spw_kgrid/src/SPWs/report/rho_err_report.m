function rho_err_report(cry,m)
fprintf('%s\n',repelem('-',20));
fprintf('Rho Error Report:\n')
rho_normal = m.rho_1;
rho_spw = m.rho_3;
delta_rho = rho_normal - rho_spw;
clear rho_spw
err_L2 = norm(delta_rho(:)) / norm(rho_normal(:));
fprintf('  Relative L2 Norm = %2.2e\n',err_L2);
err_max = max(abs(delta_rho(:)));
fprintf('  Max Absolute Error (L_Inf Norm) = %2.2e\n',err_max);
dV = cry.vol / (cry.n1 * cry.n2 *cry.n3);
err_elec = sum(abs(delta_rho(:))) * dV;
fprintf('  Integrated Absolute Difference = %2.2e\n',err_elec);
end
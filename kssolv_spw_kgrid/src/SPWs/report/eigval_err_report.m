function eigval_err_report(cry,info_store)
fprintf('%s\n',repelem('-',20));
fprintf('Total Energy Error (in Hartree):\n');
if ~isempty(info_store{2})
fprintf('Extra Bands:\n');
err_etot = abs(info_store{1}.Etot - info_store{2}.Etot);
fprintf('  Absolute Error: %2.2e\n',err_etot);
fprintf('  Absolute Error per Atom: %2.2e\n',err_etot/sum(cry.natoms));
end
fprintf('SPWs:\n');
err_etot = abs(info_store{1}.Etot - info_store{3}.Etot);
fprintf('  Absolute Error: %2.2e\n',err_etot);
fprintf('  Absolute Error per Atom: %2.2e\n',err_etot/sum(cry.natoms));





fprintf('%s\n',repelem('-',20));
fprintf('Eigval Error Report (in Hartree):\n')
nbnd0 = size(info_store{1}.Eigvals{1},1);
e00_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{1}.Eigvals,'UniformOutput',false).');
if ~isempty(info_store{2})
e10_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{2}.Eigvals,'UniformOutput',false).');
dife_c = e10_mat(:) - e00_mat(:);
fprintf('Extra Bands:\n  Maximum Absolute Error: %.2e \n',max(abs(dife_c)));
end

e11_mat = cell2mat(cellfun(@(x)x(1:nbnd0),info_store{3}.Eigvals,'UniformOutput',false).');
dife_spw = e11_mat(:) - e00_mat(:);
fprintf('SPWs:\n')
fprintf('  Maximum Absolute Error: %2.2e \n',max(abs(dife_spw)));
e_L1 = sum(abs(dife_spw))/nbnd0;
fprintf('  Mean Absolute Error(L1 Norm): %2.2e \n',max(abs(e_L1)));

e_RMSE = sqrt(sum(dife_spw.^2)/nbnd0);
fprintf('  Root Mean Square Error: %2.2e \n',max(abs(e_RMSE)));

% e_MBE = sum(dife_spw)/nbnd0;
% fprintf('  Mean Bias Error: %2.2e \n',max(abs(e_RMSE)));


end
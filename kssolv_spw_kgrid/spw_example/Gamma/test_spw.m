%  startup
% nohup matlab -nosplash -nodisplay -r "test_spw; exit" < /dev/null > sg157.log 2>&1 &
clc,clear,close all;
dateStr = datestr(now, 'yyyymmdd_HHMMSS');

% Please change the folder name.
folder_path = fullfile(pwd,'sg157-C3v');

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
rhofile_path = fullfile(folder_path,'rho_data.mat');
m = matfile(rhofile_path,'Writable',true);
tempfile = fullfile(folder_path,'temp_diary.log');
diary(tempfile);

%%
% maxNumCompThreads(48);
% nThreads = maxNumCompThreads;
% disp(['Max Threads: ', num2str(nThreads)]);

% pyenv(Version="/public/home/wangly/anaconda3/envs/pymatgen/bin/python");


cd ../../;
KSSOLV_startup; 
kssolvpptype('default')

cd spw_example;
cd structure_spg;

%---------Structure input: please select one of them.---------%

% test_sg14_natom208_Sb9S4F39_setup; 
test_sg157_natom189_Na21S7Cl_O14F3_2_setup;
% test_sg220_natom156_Cs5Ag4C8IN8_setup; 


cd ..;

cd Gamma;

energycutoff = 30;


% 4. Configure the molecule (crystal)
%  % Automatic k-points using Monkhorst-Pack grid
  autokpts = [ 1 1 1 0 0 0 ];

   cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
   'ecut', energycutoff, 'name','C', 'autokpts', autokpts, 'temperature', 0 );

   clear('natom',"i","xyzlist","C","coefs","a2","atob");
   
   fprintf('Energy Cut-off is set to %.2f Ha.\n',cry.ecut);
   fprintf('Coarse Grid Energy Cut-off is set to %.2f Ha.\n',cry.ecut_coarse);
%--------------------standard option setting------------------%
opt = setksopt;
opt.mixdim = 8;
opt.maxscfiter = 45;
opt.maxcgiter = 15;
opt.maxcgiter_c = opt.maxcgiter;
opt.force = 0;
opt.symcharge = 1;
%---------------------SPW symmetry setting--------------------%
P_fcc = 1/2*[0,1,1;1,0,1;1,1,0];
P_bcc = 1/2*[-1,1,1;1,-1,1;1,1,-1];
P_base = 1/2*[1,1,0;-1,1,0;0,0,2];
P_H = 1/3*[2,-1,-1;1,1,-2;1,1,1];
p = 1/8*[1,1,1];

%-------------Change it according to input material-----------%
S = Sym(cry); 
% If you have problem loading the pymatgen, you can load this file like:
% load('S_sg157.mat','S');

% Change the S according to the symmetry of your input structure.
% S = Sbilbaoinfo(S,14,eye(3),0*p); 
S = Sbilbaoinfo(S,157,eye(3),0*p); 
% S = Sbilbaoinfo(S,220,P_bcc,0*p); 

cry = genIBZ(cry,S);
%-------------------------------------------------------------------------%
ev_store = cell(3,1);
info_store = cell(3,1);
rho_store = cell(3,1);
init_control = 0;
nbnd0 = cry.nbnd;
%------------------options for SPW------------------%
coarse_control = [0,1,1];
spw_control = [0,0,1];
opt.test_proj = 1;

opt.eigmethod = 'davidson_qe';
opt.test_proj = 1;
opt.scftol = 1e-8;
opt.cgtol = 1e-2;

opt.scftol_coarse = 1e-3;
opt.cgtol_c = opt.cgtol;
opt_report(opt);
useprofile = 0;
opt.diag_profile = 1;
%-------------------------------------------------------------------------%
for istep = 1:3

opt.coarse = coarse_control(istep);
opt.usespw = spw_control(istep);
fprintf('==============Part%d. coarse = %d, spw = %d=============\n',istep,opt.coarse,opt.usespw);

if istep == 1
    cry.nbnd = ceil(1.00*nbnd0);
else
    cry.nbnd = ceil(1.10*nbnd0);
end

if useprofile
    profile clear;
    profile on -history;
end

[~,H,X,info] = scf_spw_2gc_initAll_1(cry,opt,S);

info_store{istep} = info;
rhoName = sprintf('rho_%d', istep);
m.(rhoName) = H.rho;

clear H X info

if useprofile
   p_data = profile('info');
   profile off; 

   html_dir = sprintf('profile_html_step%d_c%d_s%d', istep, opt.coarse, opt.usespw);
   html_ff = fullfile(folder_path,html_dir);
   mat_file = sprintf('profile_data_step%d_c%d_s%d.mat', istep, opt.coarse, opt.usespw);
   mat_ff = fullfile(folder_path,mat_file);

   profsave(p_data, html_ff);
   save(mat_ff, 'p_data');

   fprintf('>>> Profiling saved: %s\n\n', html_dir);
end

end
%------------------------------Report----------------------------------%

fprintf('==============Accuracy Report=============\n')
eigval_err_report(cry,info_store);
rho_err_report(cry,m);
if opt.diag_profile 
 fprintf('Standard SCF.\n');
 analyze_diag_profile(info_store{1}.prof_history);
 fprintf('Standard SCF with extra band.\n');
 analyze_diag_profile(info_store{2}.prof_history);
 fprintf('SPW SCF with extra band.\n');
 analyze_diag_profile(info_store{3}.prof_history);
end
clear m; delete(rhofile_path);




diary off;


dairy_fileName = fullfile(folder_path, sprintf('matlab_accuracy_%s_%s.log', opt.eigmethod, dateStr));
movefile(tempfile,dairy_fileName);


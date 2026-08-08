function [X, ev, options] = diag_david_profile(H,X,ev,options)

% ==== [探查器: 启动总外壳计时] ====
t_total_start = tic;
pd = options.prof_data;
% ==================================

% read parameters in options
tol = options.tol;
maxiter = options.maxiter;
lrot = options.lrot;
h_diag = options.h_diag;
s_diag = options.s_diag;
npol = options.npol;
% whether accurate diagnolization
btype = options.btype;

npw = size(X.psi,1)/npol;
npwx = npw;
% nvec: nbnd of X
nvec = size(X.psi,2);
% nvecx: max dimension of the reduced basis set
nvecx = nvec*2;
empty_ethr = max(5*tol,1e-5);
kdim = npw*npol;
% allocate space for iteration
psi = zeros(npwx*npol,nvecx);
hpsi = zeros(npwx*npol,nvecx);
sc = zeros(nvecx);
hc = zeros(nvecx);
vc = zeros(nvecx);
ew = zeros(1,nvecx);
conv = zeros(nvec,1);

notcnv = nvec;
nbase = nvec;
conv(:) = 0;

psi(:,1:nvec) = X.psi;

% ==== [探查器: 初始 H*X 计时与操作数统计] ====
t0 = tic;
HX = H*X;
pd.t_hpsi = pd.t_hpsi + toc(t0);
pd.n_hpsi_ops = pd.n_hpsi_ops + nvec;
% =============================================

hpsi(:,1:nvec) = 2*HX.psi;
my_n = nbase;

% ==== [探查器: 初始投影与构建小矩阵计时] ====
t0 = tic;
hc(1:nbase,1:my_n) = subtimes(psi',hpsi,nbase,my_n,kdim);
sc(1:nbase,1:my_n) = subtimes(psi',psi,nbase,my_n,kdim);
for n = 1:nbase
    hc(n,n) = complex(real(hc(n,n)),0);
    sc(n,n) = complex(real(sc(n,n)),0);
    for m = (n+1):nbase
        hc(n,m) = conj(hc(m,n));
        sc(n,m) = conj(sc(m,n));
    end
end
pd.t_proj = pd.t_proj + toc(t0);
% ============================================

if lrot
    vc(1:nbase,1:nbase) = 0;
    for n = 1:nbase
        ev(n) = real(hc(n,n));
        vc(n,n) = 1;
    end
else
    % ==== [探查器: 小矩阵特征值求解计时] ====
    t0 = tic;
    [ew(1:nvec),vc(1:nvec,1:nvec)] = diag_h(hc,sc,nbase,nvec);
    pd.t_eig = pd.t_eig + toc(t0);
    % ========================================
    ev(1:nvec) = ew(1:nvec);
end

for kter = 1:maxiter
    dav_iter = kter;
    
    % ==== [探查器: 波函数更新与预条件计时] ====
    t0 = tic;
    np = 0;
    for n = 1:nvec
        if ~conv(n)
            np = np + 1;
            if np~=n
                vc(:,np) = vc(:,n);
            end
            ew(nbase+np) = ev(n);
        end
    end
    nb1 = nbase + 1;
    my_n = nbase;
    psi(1:kdim,nb1:nb1+notcnv-1) = subtimes(psi,vc,kdim,notcnv,my_n);
    id = (1:notcnv) + nbase;
    psi(:,id) = -psi(:,id).*repmat(ew(id),npw*npol,1);
    psi(1:kdim,nb1:nb1+notcnv-1) = psi(1:kdim,nb1:nb1+notcnv-1) + subtimes(hpsi,vc,kdim,notcnv,my_n);
    
    x = repmat(2*h_diag,1,notcnv) - repmat(ew(nb1:nb1+notcnv-1),npw*npol,1).*repmat(s_diag,1,notcnv);
    denm = 0.5*(1+x+sqrt(1+(x-1).*(x-1)));
    psi(:,nb1:nb1+notcnv-1) = psi(:,nb1:nb1+notcnv-1)./denm;
    
    for n = 1:notcnv
        nbn = nbase + n;
        ew(n) = psi(:,nbn)'*psi(:,nbn);
    end
    for n = 1:notcnv
        psi(:,nbase+n) = psi(:,nbase+n)/sqrt(ew(n));
    end
    
    X_tmp = X;
    X_tmp.psi = psi(:,nb1:nb1+notcnv-1);
    pd.t_upd = pd.t_upd + toc(t0);
    % ==========================================
    
    % ==== [探查器: H*X 计时与操作数统计 (仅对未收敛向量)] ====
    t0 = tic;
    HX_tmp = H*X_tmp;
    pd.t_hpsi = pd.t_hpsi + toc(t0);
    pd.n_hpsi_ops = pd.n_hpsi_ops + notcnv;
    % =========================================================
    
    hpsi(:,nb1:nb1+notcnv-1) = 2*HX_tmp.psi;
    my_n = nbase + notcnv;
    
    % ==== [探查器: 投影与构建小矩阵计时] ====
    t0 = tic;
    hc(nb1:nb1+notcnv-1,1:my_n) = subtimes(hpsi(:,nb1:end)',psi,notcnv,my_n,kdim);
    sc(nb1:nb1+notcnv-1,1:my_n) = subtimes(psi(:,nb1:end)',psi,notcnv,my_n,kdim);
    
    nbase = nbase + notcnv;
    for n = 1:nbase
        if n>= nb1
            hc(n,n) = complex(real(hc(n,n)),0);
            sc(n,n) = complex(real(sc(n,n)),0);
        end
        for m = max(n+1,nb1):nbase
            hc(n,m) = conj(hc(m,n));
            sc(n,m) = conj(sc(m,n));
        end
    end
    pd.t_proj = pd.t_proj + toc(t0);
    % ========================================
    
    % ==== [探查器: 小矩阵特征值求解计时] ====
    t0 = tic;
    [ew(1:nvec),vc(1:nbase,1:nvec)] = diag_h(hc,sc,nbase,nvec);
    pd.t_eig = pd.t_eig + toc(t0);
    % ========================================
    
    % calculate number of the convergenced vectors
    for n = 1:nvec
        if btype(n) == 1
            conv(n) = (abs(ew(n) - ev(n)) < tol);
        else
            conv(n) = (abs(ew(n) - ev(n)) < empty_ethr);
        end
    end
    notcnv = sum(~conv);
    ev(1:nvec) = ew(1:nvec);
    
    if notcnv == 0||nbase+notcnv>nvecx||dav_iter==maxiter
        
        % ==== [探查器: 记录重启次数] ====
        if notcnv > 0  && nbase+notcnv > nvecx
            pd.n_restarts = pd.n_restarts + 1;
            pd.restart_dims(end+1) = nbase;
        end
        % ================================
        
        % ==== [探查器: 坍缩与收尾更新计时] ====
        t0 = tic;
        my_n = nbase;
        evc = subtimes(psi,vc,kdim,nvec,my_n);
        X.psi = evc;
        options.notcnv = notcnv;
        
        if notcnv == 0
            pd.t_upd = pd.t_upd + toc(t0);
            options.dav_iter = options.dav_iter + dav_iter;
            % 提前退出前的总计时结算
            pd.t_total = pd.t_total + toc(t_total_start);
            options.prof_data = pd;
            return;
        elseif dav_iter == maxiter
            pd.t_upd = pd.t_upd + toc(t0);
            options.dav_iter = options.dav_iter + dav_iter;
            % 提前退出前的总计时结算
            pd.t_total = pd.t_total + toc(t_total_start);
            options.prof_data = pd;
            return;
        end
        
        psi(:,1:nvec) = evc;
        psi(1:kdim,nvec+1:2*nvec) = subtimes(hpsi,vc,kdim,nvec,my_n);
        hpsi(:,1:nvec) = psi(:,nvec+1:2*nvec);
        nbase = nvec;
        hc(1:nbase,1:nbase) = 0;
        sc(1:nbase,1:nbase) = 0;
        vc(1:nbase,1:nbase) = 0;
        
        for n = 1:nbase
            hc(n,n) = complex(ev(n),0);
            sc(n,n) = 1;
            vc(n,n) = 1;
        end
        pd.t_upd = pd.t_upd + toc(t0);
        % ======================================
    end
end

% 正常循环结束的总计时结算
pd.t_total = pd.t_total + toc(t_total_start);
options.prof_data = pd;

end
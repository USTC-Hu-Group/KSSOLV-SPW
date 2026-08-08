function prec = genHprecSpw(H,Gprec,varargin)
% GENPREC generates the preconditioner.
%    X0 = GENPREC(H) generates the preconditioner for the Hamiltonian,
%    which is in the format of Wavefun.
%
%    BX0 = GENPREC(BH) generates the preconditioners for the Bloch
%    Hamiltonian, which is in the format of BlochWavefun.
%
%   See also Ham, BlochHam, Wavefun, BlochWavefun.
%-----------------------gpa---------------------------%
% Now the input can ONLY be BlochHam.
% This code is ONLY suitable for SINGLE k-point situation. 

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
ng = size(Gprec.dim,1);
p = cell(ng,1);
gkincell = cell(ng,1);
if isempty(varargin)
    ik = 1;
else
    ik = varargin{1};
end
idxnz = H.idxnz{ik};
for ig = 1:ng
    M = Gprec.gmat{ig,1};
    gkincell{ig} = zeros(size(M,2),1);
    if isa( H, 'BlochHam' )
        %nkpts = H.nkpts;
        %gkincell{ig} = diag(M'*diag(H.gkincell{ig})*M);
        for i = 1:size(M,2)
            gkincell{ig}(i) = sum((M(:,i) .^ 2) .* H.gkincell{ik});
        end
    else
        %nkpts = 1;
        gkincell = {H.gkin};
    end
    
    %for ik = 1:nkpts
        X  = gkincell{ig};
        Y  = 27.0 + X.*(18.0 + X.*(12.0 + 8.0*X));
        p{ig}  = Y./(Y + 16.0*X.^4);
    %end
end    
    if isa( H, 'BlochHam' )
        prec = SpwWavefunSet(p,n1,n2,n3,idxnz,H.wks);
        prec = genGpaInfo(prec,Gprec);
    else
        prec = SpwWavefun(p{1},n1,n2,n3,idxnz);
    end

end
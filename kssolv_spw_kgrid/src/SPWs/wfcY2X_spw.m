function X = wfcY2X_spw(Y)
%WFCY2X Reconstruct full-space wavefunctions from SPW coefficients.
%
% The transformation is exact:
%     X_j^(alpha) = P_j^(alpha) Y^(alpha)
%
% This implementation preallocates X.psi to avoid repeated array
% concatenation and copying.

    if ~isa(Y, 'SpwWavefunSet')
        error('wfcY2X:InvalidInput', ...
              'Input must be an SpwWavefunSet object.');
    end

    X = Wavefun();

    % -------------------------------------------------------------
    % Determine the final matrix dimensions.
    % -------------------------------------------------------------
    nrows = 0;
    ncols_total = 0;
    prototype = [];

    for i = 1:Y.nkpts
        Yi = Y{i}.gpsi;

        if isempty(Yi)
            continue;
        end

        if isempty(prototype)
            prototype = Yi;
            nrows = size(Y.gmat{i,1}, 1);
        end

        % Each irrep row produces one full-space block.
        ncols_total = ncols_total + Y.dim(i) * size(Yi, 2);
    end

    if isempty(prototype)
        X.psi = [];
    else
        % Allocate the complete output exactly once.
        X.psi = zeros(nrows, ncols_total, 'like', prototype);

        first_col = 1;

        % ---------------------------------------------------------
        % Reconstruct each partner-row block directly into its
        % final location.
        % ---------------------------------------------------------
        for i = 1:Y.nkpts
            Yi = Y{i}.gpsi;

            if isempty(Yi)
                continue;
            end

            ncols_block = size(Yi, 2);

            for j = 1:Y.dim(i)
                cols = first_col:(first_col + ncols_block - 1);

                Pij = Y.gmat{i,j};
                X.psi(:, cols) = Pij * Yi;

                first_col = first_col + ncols_block;
            end
        end
    end

    % -------------------------------------------------------------
    % Copy wavefunction metadata.
    % -------------------------------------------------------------
    X.n1        = Y{1}.n1;
    X.n2        = Y{1}.n2;
    X.n3        = Y{1}.n3;
    X.idxnz     = Y{1}.idxnz;
    X.iscompact = Y{1}.iscompact;
    X.trans     = Y{1}.trans;
    X.ispin     = Y{1}.ispin;
end
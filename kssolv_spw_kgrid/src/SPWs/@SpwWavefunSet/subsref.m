function X = subsref(BX,S)
    % BLOCHWAVEFUN/SUBSREF Subsref function for Bloch wave function class
    %    X = BX{ik} returns the a single wave function.
    %
    %    See also Wavefun, BlochWavefun.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    
    
    switch S(1).type
        case '{}'
            ik = S(1).subs{1};
            if numel(ik) > 1
                error('BlochWavefun:subsref', '仅支持单个索引值（Single index only）。')
            end
            
            % 1. 首先提取基础对象
            X = BX.spwwavefuncell{ik};
            
            % 2. 自动赋值 ik (核心修改)
            X.ik = ik;
            
            % 3. 同步赋值其他相关属性
            if ~isempty(BX.gmat)
                X.gmat = BX.gmat{ik,1};
            end

            % 如果是自旋极化计算，将全局索引 ik 转换为 k 点索引
            if BX.nspin == 2
                % 这里的逻辑假设 X.ik 已经是全局索引
                X.ik = X.ik - BX.nkpts * (X.ispin - 1);
            end
            
            X.wks = BX.wks;
            X.irep = ik;
            
            % 4. 处理链式调用 (例如 BX{ik}.psi 或 BX{ik}(m))
            if numel(S) > 1
                % 注意：此时 X 已经是赋值过 ik 的对象了
                X = builtin('subsref', X, S(2:end));
            end
            
        otherwise
            % 处理 .属性 或 () 索引
            X = builtin('subsref', BX, S);
    end
end
function X = subsref(BX,S)
 if isa(BX,'SpwprecSet')
     % disp('ff')
    switch S(1).type
        case '{}'
            ik = S(1).subs{1};
            if numel(ik) > 1
                error('Wrong sub index.')
            end
            if numel(S) > 1
                X = builtin('subsref',BX.spwpreccell{ik},S(2:end));
            else
                X = BX.spwpreccell{ik};           
            end
        otherwise
            X = builtin('subsref',BX,S);
    end
 end
end
    
function params = updateParams(params,newparams)

% update the parameters as supplied (only the supplied parameters are changed)

newfields = fieldnames(newparams);
for i = 1:length(newfields)
    if(isfield(params,newfields{i}))
        params.(newfields{i}) = newparams.(newfields{i});
    else
        %echo([newfields{i} ' is not a recognized parameter - possible typo ?????']);
    end
end

return 
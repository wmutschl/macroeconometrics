function VECMprint(VECM,opt)
    
fprintf('\n****************************')
fprintf('\nParameter Estimates\n')
        
paramNames = VECM.paramNames;
paramVals = VECM.paramVals;
fprintf('\nr = %d',opt.coint_r)
fprintf('\n------\n')

if opt.coint_r == 0 % Remove empty parameters from display
    paramNames(ismember(paramNames,{'A','B','c0','d0'})) = [];
end

for param = 1:length(paramNames)
    fprintf('%s =\n',paramNames{param})
    disp(paramVals.(paramNames{param}))
end
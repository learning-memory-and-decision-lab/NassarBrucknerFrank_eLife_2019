function [params, negLogLike]=fitLinearModWCircErrs(data)
%% GOAL: Fit a regression-like model to data on a circle (Updates) with any 
% set of predictors. 


% data.Y = ydata
% data.X = xdata
% data.includeUniform = do you want to include a uniform mixture component?


% vm concentration:


if ~isfield(data, 'startPoint')
    startPoint=[5, zeros(1, size(data.X, 2))];
    makeStartPoint=true;
else
    makeStartPoint=false;
    startPoint=data.startPoint;
end

if ~isfield(data, 'whichParams')
data.whichParams=true(size(startPoint));
end


lb=[.0001, ones(1, size(data.X, 2)).*-5];

%ub=[1000, ones(1, size(data.X, 2)).*100];
ub=[100, ones(1, size(data.X, 2)).*5];


% Set startpoint and boundaries on uniform if you include it in mix. 
if isfield(data, 'includeUniform')&data.includeUniform==1
    if makeStartPoint
    startPoint(end+1)=.5;
    end
    lb(end+1)=0;
    ub(end+1)=1;
elseif ~isfield(data, 'includeUniform')
    data.includeUniform=false;
end


%options = optimset('Algorithm','interior-point', 'MaxPCGIter', 5000, 'MaxProjCGIter', 5000, 'MaxSQPIter', 5000, 'MaxRLPIter', 5000);
options=optimset('Algorithm', 'interior-point', 'Display','final');%rasmus 26.07 %options = optimoptions('fmincon','Algorithm','interior-point');

[params] = fmincon(@circLinPredMod, startPoint(data.whichParams), [], [], [], [], lb(data.whichParams), ub(data.whichParams), [], options);


negLogLike=circLinPredMod(params);

    function [negLogLike]=circLinPredMod(params)
       %keyboard
       % warning('coeffs gerändert')
      % keyboard
       gParams(data.whichParams)=params(data.whichParams);
       gParams(~data.whichParams)=startPoint(~data.whichParams);
       params=gParams;
       
       coeffs=params(2:end-data.includeUniform);
       %coeffs = params;
       yHat=data.X*coeffs';
        allVM_likelihoods=circ_vmpdf(data.Y, yHat, params(1));
       
        % add uniform mixture component.
        if data.includeUniform==1
            allVM_likelihoods=allVM_likelihoods.*(1-params(end))+ (1./(2.*pi)).*(params(end));
        end    
        
        negLogLike=-1.*sum(log(allVM_likelihoods));
 
        % IF there is a uniform mixture, adjust the likelihoods
        % appropriately.

        if isfield(data, 'priorWidth')
            priorProb=sum(log(normpdf(coeffs, data.priorMean, data.priorWidth)));
            if priorProb <1e-300
               priorProb=1e-300; % set some minimum prior probability (flat outside some range)
            end
            negLogLike=negLogLike-priorProb;
        end
        
        if ~isfinite(negLogLike) 
            keyboard
        end

    end
end
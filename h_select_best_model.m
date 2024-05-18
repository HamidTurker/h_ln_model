%% Description
% This code will implement forward feature selection in order to determine
% the simplest model that best describes neural spiking. First, the
% highest-performing single-variable model is identified. Then, the
% highest-perfmoring double-variable model that includes the
% single-variable model is identified. This continues until the full model
% is identified. Next, statistical tests are applied to see if including
% extra variables significantly improves model performance. The first time
% that including variable does NOT signficantly improve performance, the
% procedure is stopped and the model at that point is recorded as the
% selected model.

% the model indexing scheme:
% phst, phs, pht, pst, hst, ph, ps, pt, hs, ht, st, p,  h,  s,  t
% 1      2    3    4    5    6  7   8   9   10  11  12  13  14  15

% Max number of predictors
numPreds = length(ln_params);
% Vector to save top chosen models
topModels = zeros(numPreds,1);

% Formatting of testFit's 6 columns: var ex, correlation, llh increase, mse, # of spikes, length of test data
testFit_mat = cell2mat(testFit);
LLH_values = reshape(testFit_mat(:,3),numFolds,numModels); % Grab the LLH increase values

% Find the best model
mustInclude=[];
for p = 1:numPreds
    
    % What are the p-th order models?
    % For instance, with p=1, find all models that only have one predictor
    % With p=2, find all models that have 2 predictors, and so on..
    if isempty(mustInclude) % Do we have a model selection of a lower order already?
        pModels=[];
        for m = 1:length(modelType)
            if sum(cell2mat(modelType(m)))==p
                pModels = [pModels, m];
            end
        end
        % What's the best among these p-th order models?
        [~,idx] = max(nanmean(LLH_values(:,pModels)));
        topP = pModels(idx);
        topModels(p) = topP;
        mustInclude = cell2mat(modelType(topP));

    else % We have a lower-order model already that must be included on this order
        pModels=[];
        for m = 1:length(modelType)
            if sum(cell2mat(modelType(m)))==p
                pModels = [pModels, m];
            end
        end
        % Take the candidate pModels on this order and filter to include
        % only those pModels that also include the best model of the lower
        % order, because that one must be included
        hold=[];
        for k = 1:length(pModels)
            
            % Check whether this p-th order model contains the lower order
            % model nested within it. Hold on to it, if that's the case.
            checkModel = cell2mat(modelType(pModels(k)));
            
            valid=1;
            while valid
                for l = 1:numPreds
                    if (mustInclude(l))
                        if (mustInclude(l))+checkModel(l) == 2
                            continue
                        else
                            valid=0;
                        end
                    end
                end
                break
            end
            
            if (valid)
                hold = [hold, pModels(k)];
            end

        end
        pModels = hold; % Only these models are valid
        [~,idx] = max(nanmean(LLH_values(:,pModels)));
        topP = pModels(idx);
        topModels(p) = topP;
        mustInclude = cell2mat(modelType(topP));

    end
end
% topModels now lists the best models at each order where each higher order
% has the best lower order models nested within it

% Let's compare the best fitting models against each other based on their llh values
LLH_tops = {}; % Start by finding the corresponding llh values for each model
for k = 1:numPreds
    LLH_tops{k} = LLH_values(:,topModels(k));
end

% We're gonna do a forward search (i.e., order 1 against order 2, order 2
% against order 3, order 3 against order 4, and so on. So, the number of
% comparisons equals numPreds-1.
p_llh = [];
for k = 1:(numPreds-1)
    % Compare next-higher-order against lower and save p-value
    [p_llh(k),~] = signrank(LLH_tops{k+1},LLH_tops{k},'tail','right');
end

% Move through the p-values and halt when the next model is not significant
selected_order = []; halt = 0;
while ~halt
    for k = 1:(numPreds-1)
        if (p_llh(k) > .05)
            selected_order = k; % Halt at current model order
            halt = 1;
            break
        end
        selected_order = numPreds; % Full model
    end
    halt = 1;
end
selected_model = modelType(topModels(selected_order)); % Selected model

% Re-set if selected model is not above baseline
pval_baseline = signrank(LLH_values(:,selected_order),[],'tail','right');
if pval_baseline > 0.05
    selected_model = NaN;
end





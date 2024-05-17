function [f, df, hessian] = h_ln_poisson_model(params,data,modelType,ln_params,ln_params_dims)

%% @ HBT - 2024 May 16
%    Adjusted the original code to allow for flexible number of parameters

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% Compute the firing rate
u = X * params;
rate = exp(u);

% Roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data
reg_w = 5e1;

% Start computing the Hessian
rX = bsxfun(@times,rate,X);
hessian_glm = rX'*X;

%% Find parameters and compute their roughness penalties
% Initialize parameter-relevant variables
J={}; J_g={}; J_h={};
for l=1:length(ln_params)
    J{l}=0; J_g{l}=[]; J_h{l}=[];
end

% Find the parameters for the current modelType
[params_found] = find_param(params,modelType,ln_params);

% Compute the contribution for f, df, and the hessian
for l = 1:length(params_found)
    if ~isempty(cell2mat(params_found(l)))

        if (ln_params_dims{l}) == "1d"
            [J{l},J_g{l},J_h{l}] = rough_penalty_1d(cell2mat(params_found(l)),reg_w);

        elseif (ln_params_dims{l}) == "2d"
            [J{l},J_g{l},J_h{l}] = rough_penalty_2d(cell2mat(params_found(l)),reg_w);

        elseif  (ln_params_dims{l}) == "1d circ"
            [J{l},J_g{l},J_h{l}] = rough_penalty_1d_circ(cell2mat(params_found(l)),reg_w);

        end
    end
end

%% Compute f, the gradient, and the hessian 
f = sum(rate-Y.*u) + sum(cell2mat(J));
df = real(X' * (rate - Y) + vertcat(J_g{1:end}));
hessian = hessian_glm + blkdiag(J_h{1:end});

%% Local functions called in the above script
% Find the right parameters given the model type
function [params_found] = find_param(params,modelType,ln_params)

    % Initialize
    for l=1:length(ln_params)
        params_found{l} = [];
    end
    
    % Find params
    for l=1:length(modelType)
        if l == 1
            start_idx = 1;
        else
            validsum = 0;
            for m=1:(l-1)
                if (modelType(m))
                    validsum = validsum + ln_params(m);
                end
            end
            start_idx = 1+validsum;
        end
        end_idx = start_idx+ln_params(l)-1;
        if (modelType(l))
            params_found{l} = params(start_idx:end_idx);
        end
    end
    

% Cost functions (Cost J, gradient of cost J, hessian of cost J)
function [J,J_g,J_h] = rough_penalty_2d(params,beta)

    numParam = numel(params);
    D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
    DD1 = D1'*D1;
    M1 = kron(eye(sqrt(numParam)),DD1); M2 = kron(DD1,eye(sqrt(numParam)));
    M = (M1 + M2);
    
    J = beta*0.5*params'*M*params;
    J_g = beta*M*params;
    J_h = beta*M;

function [J,J_g,J_h] = rough_penalty_1d_circ(params,beta)
    
    numParam = numel(params);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    
    % to correct the smoothing across first and last bin
    DD1(1,:) = circshift(DD1(2,:),[0 -1]);
    DD1(end,:) = circshift(DD1(end-1,:),[0 1]);
    
    J = beta*0.5*params'*DD1*params;
    J_g = beta*DD1*params;
    J_h = beta*DD1;

function [J,J_g,J_h] = rough_penalty_1d(params,beta)

    numParam = numel(params);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    J = beta*0.5*params'*DD1*params;
    J_g = beta*DD1*params;
    J_h = beta*DD1;



function [ preds ad_pred ] = bcrf_step(bcrf, func_id)

    % Determine unknown proteins
  bcrf.pred.unknown = bcrf.unknown;

    % Allocate memory
  bcrf.pred.params  = zeros(bcrf.Niter+1, 2 * bcrf.data.meas_count + 1);
  bcrf.pred.labels  = zeros(bcrf.Niter+1, bcrf.data.orf_count);
  bcrf.pred.prob1   = zeros(bcrf.Niter+1, length(bcrf.pred.unknown));
  bcrf.pred.log_r   = zeros(1, bcrf.Niter+1);
  bcrf.pred.log_u   = zeros(1, bcrf.Niter+1);
  bcrf.pred.tauc    = [];
  bcrf.pred.func_id = func_id;
  
    % To predict the true parameters
  %bcrf.pred.labels(1, :) = bcrf.data.func_assoc(:, func_id);
  %bcrf.pred.slimmed = slim(bcrf);
  %[ glm_x, glm_y ] = glmset(bcrf, bcrf.data.rel, 1);
  %bcrf.pred.true_params = glmfit(glm_x, glm_y, 'binomial','link','logit');

    % Get labels
  bcrf.pred.labels(1, :) = bcrf.data.func_assoc(:, func_id);
  bcrf.pred.labels(1, bcrf.pred.unknown) = -1;

    % Calculate slimmed network, find training
  bcrf.pred.slimmed = slim(bcrf);

  train = setdiff(bcrf.data.rel, [bcrf.pred.slimmed.H_U bcrf.pred.unknown ]);

    % Predict preliminary parameters
  bcrf.pred.params(1, :) = 1;
  if ~isempty(train)
    [ glm_x, glm_y ] = glmset(bcrf, train, 1);
    [ glm.coeff glm.dev glm.stats ] = glmfit(glm_x, glm_y, 'binomial','link','logit');
    bcrf.pred.params(1,:) = mvnrnd(glm.coeff, glm.stats.covb, 1);
  end

  occ1 = length(find(bcrf.pred.labels(1, :) == 1));
  occ0 = length(find(bcrf.pred.labels(1, :) == 0));
  bcrf.pred.p1 = 1-occ1/(occ1 + occ0);

    % Give the unknown proteins random values
  bcrf.pred.labels(1, bcrf.pred.unknown) = binornd(ones(length(bcrf.pred.unknown),1), bcrf.pred.p1);

  for i = 2:(bcrf.Niter + 1)

    [ bcrf.pred.labels(i,:) bcrf.pred.prob1(i, :) ] = bcrf_inference(bcrf, i-1);
    [ bcrf.pred.params(i,:) bcrf.pred.log_r(i) bcrf.pred.log_u(i) ] = bcrf_training(bcrf, i-1);
    fprintf('%d: %d\n', func_id, i);
  end

  bcrf.pred.tauc = zeros(2, length(bcrf.Niter));
  for i = 1:bcrf.Niter
    preds = zeros(3, length(bcrf.pred.unknown));
    preds(1, :) = bcrf.pred.unknown;
    preds(2, :) = bcrf.data.func_assoc(bcrf.pred.unknown, func_id);
    preds(3, :) = mean(bcrf.pred.prob1(i, :),1);
    
    a = auc_stair(preds);
    np        = [ i ; a ];
    bcrf.pred.tauc = [ bcrf.pred.tauc np ];
  end

  if ~isempty(bcrf.pred.unknown)
    preds = zeros(3, length(bcrf.pred.unknown));
    preds(1, :) = bcrf.pred.unknown;
    preds(2, :) = bcrf.data.func_assoc(bcrf.pred.unknown, func_id);
    preds(3, :) = mean(bcrf.pred.prob1(.2*bcrf.Niter:5:end, :),1);
  else
    preds       = zeros(3, length(bcrf.data.rel));
    preds(1, :) = bcrf.data.rel;
    preds(2, :) = bcrf.data.func_assoc(bcrf.data.rel, func_id);
    
    params      = mean(bcrf.pred.params(end-200:end,:));
    feats       = glmset(bcrf, bcrf.data.rel, 1);
    feats       = [ones(length(feats), 1) feats];
    ee          = exp(params * feats');
    preds(3, :) = ee ./ (1 + ee);
  end

  ad_pred = bcrf.pred;

end

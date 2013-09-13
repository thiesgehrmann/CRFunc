function [ params_out log_r log_u ] = bcrf_training(bcrf, iter)

  params_curr = bcrf.pred.params(iter, :);

  [ glm_x glm_y ] = glmset(bcrf, bcrf.data.rel, iter);
  [ glm.coeff glm.dev glm.stats ] = glmfit(glm_x, glm_y, 'binomial','link','logit');
  params_prop = mvnrnd(glm.coeff, glm.stats.covb, 1);
  
  feats = [ ones(1, size(glm_x, 1)) ; glm_x' ];
  
  prob1_curr = exp(params_curr * feats);
  prob1_curr = prob1_curr ./ (1 + prob1_curr);

  prob1_prop = exp(params_prop * feats);
  prob1_prop = prob1_prop ./ (1 + prob1_prop);

  glm_y = logical(glm_y);

  log_curr = sum(log(prob1_curr(glm_y)));
  log_curr = log_curr + sum(log(1 - prob1_curr(~glm_y)));
 
  log_prop = sum(log(prob1_prop(glm_y)));
  log_prop = log_prop + sum(log(1 - prob1_prop(~glm_y)));

  log_r  = min(0, log_prop - log_curr);

  log_u = log(unifrnd(0,1));
  if (log_r >= log_u) && isfinite(log_prop)
    params_out = params_prop;
  else
    params_out = params_curr;
  end

  log_r = log_curr;

end

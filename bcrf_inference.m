function [ newlab prob1 ] = bcrf_inference(bcrf, iter)

  params = bcrf.pred.params(iter, :)';
  newlab = bcrf.pred.labels(iter, :)';
  [ feats discard ] = glmset(bcrf, bcrf.pred.unknown, iter);

    % Calculate the probability of the label being 1
  feats = [ ones(1, length(bcrf.pred.unknown)) ; feats' ];
  lodd  = exp(params' * feats);
  prob1 = lodd ./ (1 + lodd);
  prob1 = max(0, prob1);

    % Sample from the distribution
  thresh = unifrnd(0, 1, 1, length(bcrf.pred.unknown));
  thresh = prob1 - thresh;
  pass   = thresh >= 0;

    % assign labels
  unlab       = zeros(length(bcrf.pred.unknown), 1);
  unlab(pass) = 1;
  newlab(bcrf.pred.unknown) = unlab;
  
  prob1(~isfinite(prob1)) = 0;

end

function [ glm_x glm_y ] = glmset(bcrf, prots, iter)

  slim_complete = bcrf.pred.slimmed.prec;
  H_U           = bcrf.pred.slimmed.H_U;
  N_U           = bcrf.pred.slimmed.N_U;
  N_E           = bcrf.pred.slimmed.N_E;
  labels = bcrf.pred.labels(iter,:)';

    % For each relevant protein
  for prot = intersect(H_U, prots)

      % If it is, then sum over the measurements given the label
    neighlabs = bcrf.pred.labels(iter, N_U{prot});
    neighlabs = logical(neighlabs);
    neighev   = N_E{prot};

    neighpos = sum(neighev(neighlabs,:), 1);
    neighneg = sum(neighev(~neighlabs,:), 1);

    slim_complete(prot, :) = slim_complete(prot, :) + [ neighpos neighneg ];

  end

    % Return glm_y = f(glm_x)
  glm_y = labels(prots);
  glm_x = slim_complete(prots, :);

end

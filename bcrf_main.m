function [ preds ] = bcrf_main(data, rel, unknowns, iter, func)

      % Select subnetwork
  if ~isempty(rel)
    data.rel = rel;
  end

  bcrf.data    = data;
  bcrf.unknown = unknowns;
  bcrf.Niter   = iter;

  preds = bcrf_step(bcrf, func);

end


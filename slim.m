function [ slimmed ] =  slim(bcrf)

  fprintf('Performing slimming\n');
  
  labels          = bcrf.pred.labels(1, :)';
  slimmed.prec    = zeros(bcrf.data.orf_count, 2*bcrf.data.meas_count);
  slimmed.neighev = cell(bcrf.data.orf_count, 1);
  slimmed.N_U     = cell(bcrf.data.orf_count, 1);
  slimmed.H_U     = [];
  slimmed.N_E     = cell(bcrf.data.orf_count, 1);
  
  for prot = bcrf.data.rel

    neighev = permute(bcrf.data.evidence(:, prot, :), [ 1 3 2 ]);
    empty   = sum(neighev, 2) == 0;

      % Get a list of neighboring unknowns
    neighun = find(labels .* ~empty == -1)';

      % Sum over the known neighbors for positive and negative features
    empty(bcrf.pred.unknown)      = 1;
    neighev(bcrf.pred.unknown, :) = 0;

      % Select just the non empty rows
    neighev = neighev(~empty,:);
    rellabs = logical(labels(~empty));

      % Sum over the rows for each label
    neighpos = sum([ neighev(rellabs, :)  ; zeros(1, bcrf.data.meas_count) ], 1);
    neighneg = sum([ neighev(~rellabs, :) ; zeros(1, bcrf.data.meas_count) ], 1);

      % Return

    slimmed.prec(prot, :) = [ neighpos neighneg ];
    if ~isempty(neighun)
        % Proteins we have to consider
      slimmed.H_U       = [ slimmed.H_U prot ];
        % The evidence we need
      slimmed.N_E{prot} = permute(bcrf.data.evidence(prot, neighun, :), [ 2 3 1 ]);
        % The unknown neighbors of the protein
      slimmed.N_U{prot} = neighun;
    end

  end
 
end

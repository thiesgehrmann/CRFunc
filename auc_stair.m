function [ auc fpr tpr ] = auc_stair(preds)

  prot_lab = preds(2, :);
  prot_pr1 = preds(3, :);

  pos_id = find(prot_lab == 1);
  neg_id = find(prot_lab == 0);

  tpr = zeros(1, length(pos_id));
  fpr = zeros(1, length(pos_id));

  pos = sort(prot_pr1(pos_id), 'descend');
  neg = sort(prot_pr1(neg_id), 'descend');
  
  for prot_i = 1:length(pos_id)
    tpr(prot_i) = prot_i;
    fpr(prot_i) = length(find(neg >= pos(prot_i)));
  end

  tpr = tpr / length(pos_id);
  fpr = fpr / length(neg_id);

  tpr = [ 0 tpr 1 ];
  fpr = [ 0 fpr 1 ];

  tpr = reshape(repmat(tpr, 2,1), 1, length(tpr)*2);
  fpr = reshape(repmat(fpr, 2,1), 1, length(fpr)*2);

  tpr = tpr(1:end-1);
  fpr = fpr(2:end);

  auc = trapz(fpr, tpr);

end

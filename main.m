load('data/hulsmanet_PPI_GI_CC_union.mat');

bcrf.data = data;
bcrf.Niter = 200;
reps = 1;

func_id = 973;

mask_pos = ceil(.2*length(find(data.func_assoc(data.rel,func_id) == 1))) + 2;
mask_neg = 300 - mask_pos;

auc = zeros(reps);


for i = 1:reps
  
  bcrf.unknown = [ randsample(data.rel(data.func_assoc(data.rel, func_id) == 1), mask_pos) ...
                   randsample(data.rel(data.func_assoc(data.rel, func_id) == 0), mask_neg) ];
  
  [ preds pred_ad ] = bcrf_step(bcrf, func_id);
  auc(i)  = auc_stair(preds);

end

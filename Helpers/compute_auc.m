function auc = compute_auc(probs,recVec)
% computes AUC. I had been using perfcurc(), but the rest of the
% team does it this way with log spacing. I don't understand why,
% but doing it this way for consistency.

hr = NaN(1,100);
fa = NaN(1,100);
vals = logspace(log10(.000001),log10(1),100);
for i = 1:length(vals)
   hr(i) = sum(probs(recVec==1)>vals(i))/sum(recVec==1);
   fa(i) = sum(probs(recVec==0)>vals(i))/sum(recVec==0);
end
auc = trapz(fliplr(fa),fliplr(hr));
function error = compute_error(h,w1,w2)
if length(w1) == length(w2)
%2-norm
    error = sqrt(h)*norm(w1-w2,2);%sqrt(sum((w1-w2).^2 ));
%1-norm
%error = norm(w1-w2,Inf);
%error = h*norm(w1-w2,1);%sqrt(sum((w1-w2).^2 ));
%    error = sqrt(h)*norm(w1-w2,2);%sqrt(sum((w1-w2).^2 ));
else
%    wt = mean(reshape(w2,2,length(w2)/2))';
    it = 1 : 2 :length(w2)-1;
    it1 = it+1;
    wt = .25*(w2(it,it)+w2(it1,it)+w2(it,it1)+ w2(it1,it1));

    error = sqrt(h)*norm(w1-wt,2);
end
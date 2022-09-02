function DD=initial(p,N)
lb = ones(1,N); % lower bounds for A,V,h,l and b
ub = (p+1)*ones(1,N); % upper bounds for A,V,h,l and b
X = lhsdesign(p,N,'criterion','correlation');
D = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)))-0.5;
DD = D';
end
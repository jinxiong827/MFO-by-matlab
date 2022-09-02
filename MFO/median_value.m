function a=median_value(A)
for i=1:length(A)
    B(i)=sum(abs(log(A/A(i))));
end
[~,b]=min(B);
a=A(b);
end
function d=distance(A,B)
if length(A)~=length(B)
    error('需要保证输入数列长度相等')
end
d=0;
for i=1:length(A)
   d=d+(A(i)-B(i))^2;
end
d=sqrt(d);   
end
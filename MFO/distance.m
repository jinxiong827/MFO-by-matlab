function d=distance(A,B)
if length(A)~=length(B)
    error('��Ҫ��֤�������г������')
end
d=0;
for i=1:length(A)
   d=d+(A(i)-B(i))^2;
end
d=sqrt(d);   
end
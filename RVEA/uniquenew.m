function [NewDec2,index]=uniquenew(NewDec,A0dec,d)
k=1;
if size(NewDec,1)~=0
for i=1:size(NewDec,1)
    err=min(sum((A0dec-repmat(NewDec(i,:),size(A0dec,1),1)).^2,2).^0.5);
    if err>d
        NewDec2(k,:)=NewDec(i,:);
        index(k)=i;
        k=k+1;
    end
end
else
    NewDec2=[];
    index=0;
end
end
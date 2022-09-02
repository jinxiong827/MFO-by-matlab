function AA=arcpoint(A,B,num,degree1)
    degree2=(num+1)*degree1;
    R=((A(1)-B(1))^2+(A(2)-B(2))^2+(A(3)-B(3))^2)^0.5/2;
    C=(A+B)/2;
%     plot3(A(1),A(2),A(3),'y*')
%     hold on
%     plot3(B(1),B(2),B(3),'y*')
%     hold on    
%     plot3(C(1),C(2),C(3),'y*')
%     hold on
    A=A-C;
    B=B-C;
    Aa(1)=(A(1)*cosd(-90)-A(3)*sind(-90))/tand(degree2/2)+C(1);
    Aa(2)=A(2)/tand(degree2)+C(2);
    Aa(3)=(A(3)*cosd(-90)+A(1)*sind(-90))/tand(degree2/2)+C(3);
%     plot(Aa(1),Aa(3),'b*')
%     hold on
    A=A+C-Aa;
    B=B+C-Aa;
    for i=1:num
        AA(1,i)=A(1)*cosd(degree1*i)-A(3)*sind(degree1*i)+Aa(1);
        AA(2,i)=A(2)+Aa(2);
        AA(3,i)=A(3)*cosd(degree1*i)+A(1)*sind(degree1*i)+Aa(3);
%         plot3(AA(1,i),AA(2,i),AA(3,i),'r*')
%         hold on
    end  
    
   
end
function angle=included_angle(A1,A2,B1,B2)
    A0=A2-A1;
    B0=B2-B1;
    angle=acosd((A0(1)*B0(1)+A0(3)*B0(3))/(A0(1)^2+A0(3)^2)^0.5/(B0(1)^2+B0(3)^2)^0.5);
end
function blade_xy=blade_profile(meanline_xy,b,max_t,profile)
max_t=max_t*b;
xtb=linspace(0,b,size(meanline_xy,2));
xt1=linspace(0,1,size(meanline_xy,2));
for i=1:length(xt1)
    switch profile
        case 'NACA'
            yt(i)=(max_t/0.2)*(0.2969*xt1(i)^0.5-0.1260*xt1(i)-0.3516*xt1(i)^2+0.2843*xt1(i)^3-0.1015*xt1(i)^4);% NACA00t profile
        case'C4'
            yt(i)=(max_t/0.2)*(0.3048*xt1(i)^0.5-0.0914*xt1(i)-0.8614*xt1(i)^2+2.1236*xt1(i)^3-2.9163*xt1(i)^4+1.9744*xt1(i)^5-0.5231*xt1(i)^6);% C4 profile
        otherwise
            error('no data')
    end
end
for i=2:size(meanline_xy,2)
    meanline_theta(i)=atand((meanline_xy(2,i)-meanline_xy(2,i-1))/(meanline_xy(1,i)-meanline_xy(1,i-1)));
end
meanline_theta(1)=-90;
for i=1:size(meanline_xy,2)
    blade_profile_xu(i)=meanline_xy(1,i)+yt(i)*cosd(meanline_theta(i)+90);        % upper surface
    blade_profile_yu(i)=meanline_xy(2,i)+yt(i)*sind(meanline_theta(i)+90);
    blade_profile_xl(i)=meanline_xy(1,i)+yt(i)*cosd(meanline_theta(i)-90);        % lower surface
    blade_profile_yl(i)=meanline_xy(2,i)+yt(i)*sind(meanline_theta(i)-90);
    if norm([blade_profile_xu(i) blade_profile_yu(i)])<norm([blade_profile_xl(i) blade_profile_yl(i)])
        [blade_profile_xl(i),blade_profile_xu(i)]=deal(blade_profile_xu(i),blade_profile_xl(i));
        [blade_profile_yl(i),blade_profile_yu(i)]=deal(blade_profile_yu(i),blade_profile_yl(i));
    end
    
end

% plot(meanline_xy(1,:),meanline_xy(2,:),'r.')
% hold on
% plot(blade_profile_xu,blade_profile_yu,'g.')
% hold on
% plot(blade_profile_xl,blade_profile_yl,'b.')
% for i=1:size(meanline_xy,2)
%     hold on
%     plot([blade_profile_xu(i) blade_profile_xl(i)],[blade_profile_yu(i) blade_profile_yl(i)])
% end

blade_xy=[blade_profile_xu;blade_profile_yu;blade_profile_xl;blade_profile_yl];
end
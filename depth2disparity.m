% % 
% % Pedro Pinto: 2015
% % updated on: 30/October
% % 
% %------------------------------------------------------------------------
function Disp_Laser = depth2disparity(T,Ima3D)

%   Disparity Image (uint16) from Laser Depth Map

baseline=-T.P1(1,4)/T.P1(1,1); %baseline
fLength=T.P2(1,1); %focal length
%Z(depth) = (focalLength * baseline) / disparity
%disparity=(focalLength * baseline)/depth
[m,n] = size(Ima3D);
disp_est = zeros(m,n);

for u=1:m
    for v=1:n   
        disp_est(u,v) = ((fLength)*baseline)/Ima3D(u,v);    
    end
end

D = double(disp_est);

I = D*256;
I(I==Inf)=-1;
I(D==0) = 1;
I(I<0) = 0;
I(I>65535) = 0;
Disp_Laser = uint16(I);

end


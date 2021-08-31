function [A,B]=runEvaluationForOne(id,grid_,method_,threshold,epsthresh)
%clear variables; close all;
 

 
data_path = 'F:\Dataset_Kitti\training\'; %Path to kitti dataset
base_dir = [data_path,'velodyne'];
calib_dir = [data_path,'calib'];
path_ima = [data_path,'image_2'];
 
%pcd_path='results/pcd/';
depth_path='results/disp_0';
%s3='results/review/';
 
%cd './code';
%mex('applyMethod.cpp')
%cd ..
 
 
 
 
%ISR-AIWALKER::HMI -4.4573e-17;0.058299
%numFramesi=0;
%numFramese=99;
 
method=method_;
grid=grid_;
%Methods List:
%min: method==0
%mean: method==1
%median: method==2
%IDW: method==3
%BF: method==4
%HMF: method==5
%HMF: method==6
%HMF: method==7
%-----Incomplete ==8
%HMF: method==9
% Bilateral==11
 
i=id;
%for i=numFramesi:1:numFramese
     
     
    %Load data
     
     
%     fd = fopen( [path_ima,ima(i).name] );
%     if fd < 1
%         fprintf('Cound not open RGB image !!!\n');    keyboard
%     else
        ImaRGB = imread( sprintf('%s/%06d.png',path_ima,i) );
%     end
%     fclose(fd);
     
    T=loadCalibration(sprintf('%06d.txt',i),calib_dir);
     
    fd = fopen(sprintf('%s/%06d.bin',base_dir,i),'rb'); %OBJECT
    if fd < 1
        fprintf('No LIDAR files !!!\n');
        keyboard
    else
        velo = fread(fd,[4 inf],'single')';
        fclose(fd);
    end
     
    % remove all points behind image plane (approximation)
    idx = velo(:,1)<5;
    velo(idx,:) = [];
     
    % project to image plane (exclude luminance)
    px = (T.P2 * T.R0_rect * T.Tr_velo_to_cam * velo')';
    px(:,1) = px(:,1)./px(:,3);
    px(:,2) = px(:,2)./px(:,3);
    % % -----------------------------------------------------------------------
    ix = px(:,1)<0.5;                 px(ix,:)=[];
    ix = px(:,2)<0.5;                 px(ix,:)=[];
    ix = px(:,1)>size(ImaRGB,2)+0.5;    px(ix,:)=[];
    ix = px(:,2)>size(ImaRGB,1)+0.5;    px(ix,:)=[];
    % % Ordering
%     Pts = zeros(size(px,1),4);
    Pts = sortrows(px,2);
     
    % % ======================= Interpolation / Upsampling :::
    c_px = floor(min(px(:,2)));
    i_size = size(ImaRGB(c_px:end,:,1));
    Ima3D = zeros( size(ImaRGB(:,:,1)) );
     
    % Simply type: mex fun_dense3D.cpp
    tic;
    Ima3D(c_px:end,:) = applyMethodCP(Pts,[c_px i_size method grid threshold epsthresh double(i)]); % MEX-file
    toc;
    % % -----------------------------------------------------------------------
    % Normalization 8 bits
    %Ima_Range = uint8(255*Ima3D/max(max(Ima3D))); % :)
      Ima_Range=Ima3D;
 
      
     
    Disp_Laser =  depth2disparity(T,Ima_Range);
     
     
A=Ima3D;
B=Disp_Laser;
 
 
 
end
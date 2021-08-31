clear all
close all
clc
st.x_min = 5;   % movement direction
st.x_max = 80;
st.y_min = -80; % right
st.y_max = 80;  % left
st.bias = 1.73; % velodyne elevation 

for p=1:7481
    frame=p-1;
    [A,B]=runEvaluationForOne(frame,13,4,0.15,0.01);
    
    dm = (st.x_max*(A-st.x_min))./(A*(st.x_max-st.x_min)); 
    dm(dm < 0) = 0; dm(dm > 1) = 1;
    dm_final  = uint8(255*dm);
  
    %% Salvando as imagens
    I=imshow(dm_final);
    Idata=I.CData;
    nome=sprintf('%s%06d','',frame);
    salvo_como='.png';
    juntar=strcat(nome,salvo_como);
    imwrite(Idata,juntar,'png');
    close all
end
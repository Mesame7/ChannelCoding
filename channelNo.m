%comments
%error 8-9 >> err4/5 >>>> err89+[7eta zyada]
%throughput
%

close all
clear all
clc
%perr=0.001

counterSentData=0;
counterNeededData=0;
perr=0.001
obj=VideoReader('foreman.avi');%video
a=read(obj);%bytes
frames=get(obj,'NumberOfFrames');
for i=1:frames
    I(i).cdata=a(:,:,:,i); % getting ith frame's bytes
end
%BefEncod is cdata
s=size(I(1).cdata);
mov(1:frames) =struct('cdata', zeros(s(1),s(2), 3, 'uint8'),'colormap', []);

for i=1:frames
    %Red Components of the Frame
    R=I(i).cdata(:,:,1);
    %Green Components of all Frames
    G=I(i).cdata(:,:,2);
    %Blue Components of the Frame
    B=I(i).cdata(:,:,3);
    Rdouble = double(R);
    Gdouble = double(G);
    Bdouble = double(B);
    Rbinx(i).bindata = de2bi(Rdouble);
    Rbin(i).bindata=reshape(Rbinx(i).bindata,1,202752);
    
    
    Gbinx(i).bindata = de2bi(Gdouble);
    Gbin(i).bindata=reshape(Gbinx(i).bindata,1,202752);
    
    
    Bbinx(i).bindata = de2bi(Bdouble);
    Bbin(i).bindata=reshape(Bbinx(i).bindata,1,202752);
end

for i=1:frames
    BefEncod(i).Fdata=[Rbin(i).bindata,Gbin(i).bindata,Bbin(i).bindata];
    
end
for i=1:frames
    BefEncod(i).Fdata1024=reshape(BefEncod(i).Fdata,[594 1024]);
end
trellis = poly2trellis(7,[171,133 ]);
PuncturingRule=[[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];[1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0];[1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0];  [1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0];[1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0]];


for j=1:30
    framePacket=BefEncod(j).Fdata1024;
    
    for i=1:594
        
        Decoded(i,1:1024)=bsc( framePacket(i,:),perr);
        counterSentData=counterSentData+1024;
        counterNeededData=counterNeededData+1024;        
        
        
    end
    dec(j,1:594,1:1024)=Decoded;
    % Rec(j).rec=A;
end



%FramesMixedColorsRecieved(m).data=
%reshape(sometempArr?,1,?608256?);
for i=1:30
    A=dec(i,:,:);
    tempAA=reshape(A,1,608256);
    SomeArrayForTesting(i,:)= tempAA;
    temp=SomeArrayForTesting(i,:);
    Redbits=temp(1:202752);
    Greenbits=temp(202753:405504);
    Bluebits=temp(405505:608256);
    RedtempCol=reshape(Redbits,25344,8);
    GreentempCol=reshape(Greenbits,25344,8);
    BluetempCol=reshape(Bluebits,25344,8);
    RedDoub=bi2de(RedtempCol);
    GreenDoub=bi2de(GreentempCol);
    BlueDoub=bi2de(BluetempCol);
    Red8= uint8(255 * mat2gray(RedDoub));
    Green8=uint8(255 * mat2gray(GreenDoub));
    Blue8=uint8(255 * mat2gray(BlueDoub));
    REDTEMPP=reshape(Red8,144,176);
    Out(i).Red=REDTEMPP;
    GREENTEMPP=reshape(Green8,144,176);
    Out(i).Green=GREENTEMPP;
    BLUETEMPP=reshape(Blue8,144,176);
    
    Out(i).Blue=BLUETEMPP;
end
for i=1:30
   mov(1,i).cdata(:,:,1) =  Out(i).Red;
    mov(1,i).cdata(:,:,2) =  Out(i).Green;
    mov(1,i).cdata(:,:,3) =  Out(i).Blue;
 
    
end
Throu=counterNeededData./counterSentData
%v = VideoWriter('NewVideo.avi',);
v = VideoWriter('nocoding_Err0001.avi','Motion JPEG AVI');
open(v);
writeVideo(v,mov);
close(v)
implay('nocoding_Err0001.avi')

%
%tempfor1024=[]
%for i =1:frames
%    for j=1:594 %144*176*3*8/1024
%   tempfor1024=cat(1,tempfor1024,BefEncod(i).Fdata(((j-1)*1024+1):j*1024)) ;
%   disp([i j])
%    end
%    BefEncod(i).F1024=tempfor1024;
%end

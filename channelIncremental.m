close all
clear all
clc

perr=0.00001:0.02:0.20001 % to show the graphs
%perr=[0.1]; % to show the video

obj=VideoReader('foreman.avi');%video
a=read(obj);%bytes
% To get the number of frames in the video you can use the following:
frames=get(obj,'NumberOfFrames');
% To extract the frames of the video so you can work on them:
for i=1:frames
    I(i).cdata=a(:,:,:,i); % getting ith frame's bytes
end
%BefEncod is cdata
% In this code, you have to generate a new video with the same size as
% the original video so you can add these line :
s=size(I(1).cdata);
mov(1:frames) =struct('cdata', zeros(s(1),s(2), 3, 'uint8'),'colormap', []);

%we now have the frames of the video consisting of red blue green
%seperate them using the next loop



% this is to sepearte each colour in a frame in 3 arrays

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
    Rbin(i).bindata=reshape(Rbinx(i).bindata,1,202752);%final 1D array of frames red colour
    
    
    Gbinx(i).bindata = de2bi(Gdouble);
    Gbin(i).bindata=reshape(Gbinx(i).bindata,1,202752);%final 1D array of frames green colour
    
    
    Bbinx(i).bindata = de2bi(Bdouble);
    Bbin(i).bindata=reshape(Bbinx(i).bindata,1,202752); %final 1D array of frames Blue colour
end
% loop to have 1D array containing the 3 colours in sequence
%(red,green,blue)
for i=1:frames
    BefEncod(i).Fdata=[Rbin(i).bindata,Gbin(i).bindata,Bbin(i).bindata];
    
end
%loop for reshaping the frames data(Fdata) inside the beforeEncoding array into a
%1024 columns , 594 rows array to occupy the (608256 bit data of the RGB
%bit stream) inside Fdata1024 attribute of the array
for i=1:frames
    BefEncod(i).Fdata1024=reshape(BefEncod(i).Fdata,[594 1024]);
end

%we now have the video bit stream as an array of RGB colours in sequence ,
%so next step to apply convolutional channel coding scheme
trellis = poly2trellis(7,[171,133 ]);
%each index in the puncRule array contain different punct. matrices for the
%incremntal convolution code
PuncturingRule=[[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];[1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0];[1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0];  [1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0];[1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0]];
counterSentData=zeros(1,max(size(perr)));
counterNeededData=zeros(1,max(size(perr)));
counterErrored=zeros(1,max(size(perr)));


for k=1:max(size(perr))
    disp(perr(k))
    for j=1:30
        disp(counterErrored(k));
        
        framePacket=BefEncod(j).Fdata1024;
        % we take each 1024 bit packet and operate on it
        for i=1:594
            
            %disp([j i]);
            
            %encode the information bits using 1/2 conv encoder then apply punc
            %rule on it to have much higher data rate= 8/9
            %and continue sending with this high rate till we find errors
            
            A(i).encoded12= convenc ( framePacket(i,:), trellis,PuncturingRule(1,:));
            A(i).encoded89=convenc ( framePacket(i,:), trellis,PuncturingRule(5,:));
            %send encoded message on BSC channel with prob error= perr
            A(i).errored89=bsc( A(i).encoded89,perr(k));
            %decode the message(with rate 8/9) after it passed by the BSC
            decoded= vitdec (A(i).errored89 ,trellis,35 ,'trunc','hard',PuncturingRule(5,:));
            
            %------ %8/9 --> 4/5-------
            if(framePacket(i,:)==decoded)%check if the message was sent without errors with the current rate
                
                %    disp('89');
                Decoded(i,1:1024)=decoded;
                
                counterSentData(k)=counterSentData(k)+1152;%used for throuput calcualtion (original+redundant message bits)
                counterNeededData(k)=counterNeededData(k)+1024; %used for throuput calcualtion (original message bits)
                errors=xor(Decoded(i,1:1024),framePacket(i,:));%calculate the number error for each packet
                errs=sum(errors(:) == 1); 
                counterErrored(k)=counterErrored(k)+errs;
                %A(i).decoded=decoded;
                continue;
                
            else %if an error occured try the next rate (4/5)
                A(i).encoded45=convenc ( framePacket(i,:), trellis,PuncturingRule(4,:));
                A(i).errored45=bsc( A(i).encoded45,perr(k));
                decoded= vitdec (A(i).errored45 ,trellis,35 ,'trunc','hard',PuncturingRule(4,:));
                if(framePacket(i,:)==decoded)%check if the message was sent without errors with the current rate
                    %4/5 --> 2/3
                    %    disp('45');
                    Decoded(i,1:1024)=decoded;
                    errors=xor(Decoded(i,1:1024),framePacket(i,:));
                    errs=sum(errors(:) == 1);
                    counterErrored(k)=counterErrored(k)+errs;
                    counterSentData(k)=counterSentData(k)+1280;  %used for throuput calcualtion (original+redundant message bits)
                    counterNeededData(k)=counterNeededData(k)+1024;
                    
                    %A(i).decoded=decoded;
                    continue;
                    
                else %if an error occured try the next rate (2/3)
                    A(i).encoded23=convenc ( framePacket(i,:), trellis,PuncturingRule(3,:));
                    A(i).errored23=bsc( A(i).encoded23,perr(k));
                    decoded= vitdec (A(i).errored23 ,trellis,35 ,'trunc','hard',PuncturingRule(3,:));
                    %2/3 --> 4/7
                    if(framePacket(i,:)==decoded) %check if the message was sent without errors with the current rate
                        
                        %disp('23');
                        Decoded(i,1:1024)= decoded;
                        errors=xor(Decoded(i,1:1024),framePacket(i,:));
                        errs=sum(errors(:) == 1);
                        counterErrored(k)=counterErrored(k)+errs;
                        counterSentData(k)=counterSentData(k)+1536; %used for throuput calcualtion (original+redundant message bits)
                        counterNeededData(k)=counterNeededData(k)+1024;
                        
                        %A(i).decoded=decoded;
                        continue;
                        
                    else%if an error occured try the next rate (4/7)
                        A(i).encoded47=convenc ( framePacket(i,:), trellis,PuncturingRule(2,:));
                        A(i).errored47=bsc( A(i).encoded47,perr(k));
                        decoded= vitdec (A(i).errored47 ,trellis,35 ,'trunc','hard',PuncturingRule(2,:));
                        %4/7 --> 1/2
                        if(framePacket(i,:)==decoded) %check if the message was sent without errors with the current rate
                            
                            % disp('47');
                            counterSentData(k)=counterSentData(k)+1792; %used for throuput calcualtion (original+redundant message bits)
                            counterNeededData(k)=counterNeededData(k)+1024;
                            
                            Decoded(i,1:1024)= decoded;
                            errors=xor(Decoded(i,1:1024),framePacket(i,:));
                            errs=sum(errors(:) == 1);
                            counterErrored(k)=counterErrored(k)+errs;
                            %A(i).decoded=decoded;
                            continue;
                            
                        else
                            
                            %disp('12');
                            %if an error occured try the next rate (1/2) aka: no puncturing
                            A(i).errored12=bsc( A(i).encoded12,perr(k));
                            decoded=vitdec (A(i).errored12 ,trellis,35 ,'trunc','hard',PuncturingRule(1,:));
                            %A(i).decoded= decoded;
                            Decoded(i,1:1024)=decoded;
                            counterSentData(k)=counterSentData(k)+2048;
                            counterNeededData(k)=counterNeededData(k)+1024;
                            errors=xor(Decoded(i,1:1024),framePacket(i,:));
                            errs=sum(errors(:) == 1);
                            counterErrored(k)=counterErrored(k)+errs;
                        end
                    end
                    
                end
            end
            
            
        end
        dec(j,1:594,1:1024)=Decoded;
        % Rec(j).rec=A;
    end
    
end
codedBitProb=counterErrored./counterNeededData
%the final value of the throuput of the system
Throu=counterSentData./counterNeededData
figure();
plot(perr,codedBitProb);
title('Error Incremental');
hold on
figure();
plot(perr,Throu);
title('Throughput Incremental');



%FramesMixedColorsRecieved(m).data=
%reshape(sometempArr?,1,?608256?);
%we have the frames after channel coding, we want to reconstruct the video
%again
for i=1:30
    %we here extract each frame then extract each colour of the RGB
    A=dec(i,:,:);
    tempAA=reshape(A,1,608256);
    SomeArrayForTesting(i,:)= tempAA;
    temp=SomeArrayForTesting(i,:);
    
    %divide the 608256 bits into 3 seperate groups again
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
%this loop is to collect the 3 RGB values after reshaping and re assembling
%a new video
for i=1:30
    mov(1,i).cdata(:,:,1) =  Out(i).Red;
    mov(1,i).cdata(:,:,2) =  Out(i).Green;
    mov(1,i).cdata(:,:,3) =  Out(i).Blue;
    
    
end
%v = VideoWriter('incErr01.avi','Uncompressed AVI');
%open(v);
%writeVideo(v,mov);

%close(v);
%implay('incErr01.avi');

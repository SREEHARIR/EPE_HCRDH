clc;
close all;
clear;

% fprintf('program Start');
% 
% Ke = 12345678;
% Kw = 12345;

% data = ['This is the message to be embedded. '...
%     'This is the continuation of the same message.\n'...
%     'This is a really long message to check the amount of problems in the word hiding. \n'...
%     'The size of the message is limited to around 512*512/8 characters in the string.\n'...
%     'This string should be fully decipherable and should be shown in the command Window.\n\n'];
% data = [data data data data data data data data data data data data data data data];

% img = imread('./images/hands1.jpg');
% img = imread('./images/satellite.png');
% img = imread('./images/airplane.pgm');

[file,path]=uigetfile('*.pgm','Select image file');
ss=strcat(path,file);
img=imread(ss);

%% Image Standardisation %%

% I = rgb2gray(img);
I = img;
I = imresize(I,[512 512]);
I = double(I);
figure('Name','Original Image','NumberTitle','off');
imshow(uint8(I));

%% Preprocess Prediction Error %%

Iinv = bitset(I,8,~bitget(I,8));
% figure('Name','I inverse');
% imshow(uint8(Iinv));

In=zeros(512);
error = zeros(512);

% z=1;
% eqpre(:,1) = [0;0];
for ii = 1:512
    for jj = 1:512
%         if ii == 1 && jj <= 16
%             continue;
        if ii == 1 && jj == 1
            delta = 0;
            deltainv = 1;
        elseif ii == 1 && jj ~= 1
            delta = abs(I(ii,jj)-In(ii,jj-1));
            deltainv = abs(Iinv(ii,jj)-In(ii,jj-1));
        elseif ii ~= 1 && jj == 1
            delta = abs(I(ii,jj)-In(ii-1,jj));
            deltainv = abs(Iinv(ii,jj)-In(ii-1,jj));
        else
            delta = abs(I(ii,jj) - floor((In(ii-1,jj) + In(ii,jj-1))/2));
            deltainv = abs(Iinv(ii,jj)- floor((In(ii,jj-1)+ In(ii-1,jj))/2));
        end
        if(delta > deltainv)
%             In(ii,jj) = Iinv(ii,jj);
            In(ii,jj) = I(ii,jj);
            error(ii,jj) = 1;
            if ii == 1 && jj <= 16
%                 In(ii,jj) = I(ii,jj) + 64*(((I(ii,jj)<Iinv(ii,jj))*2)-1);
                In(ii,jj) = I(ii,jj) + ((1+delta-deltainv)*(((I(ii,jj)<Iinv(ii,jj))*2)-1));
                error(ii,jj) = 0;
            end
            if mod(jj,8) == 0 && all(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])
                In(ii,jj) = I(ii,jj) + ((1+delta-deltainv)*(((I(ii,jj)<Iinv(ii,jj))*2)-1));
                error(ii,jj) = 0;
            end
        elseif delta == deltainv
            diff = (I(ii,jj) - Iinv(ii,jj))/abs(I(ii,jj) - Iinv(ii,jj));
            In(ii,jj) = I(ii,jj) - diff;
%             fprintf('[%i %i %i]\n',ii,jj,diff);
%             eqpre(:,z) = [ii jj];
%             z = z+1;
        else
            In(ii,jj) = I(ii,jj);
        end
    end
end

% figure('Name','Prediction Corrected vs Original');
% imshow(In~=I,[]);
% figure('Name','Error');
% imshow((error),[]);
%% Encryption %%


Ke = inputdlg({'Enter the Key for Encrypting Image'},'Image Encryption Key',[1,35],{'123456789'});
Ke = str2double((cell2mat(Ke)));
seed = Ke;
rng(seed,'twister');
S = randi(255,512);

Ie = bitxor(S,double(uint8(In)));
% figure('Name','Encrypted Image','NumberTitle','off');
% imshow(uint8(Ie));

%% Create Flag bits to surround error bits

flag = zeros([512 512]);
for ii = 1:512
    jj = 1;
    while jj <= 512
        if ii == 1 && jj <= 16
            jj = jj + 1;
            continue;
        end
        if error(ii,jj) == 1
            if flag(ii,jj) == 0
                % flag prev = 1
                if floor((jj-1)/8) ~= 0 %% current_block ~= 0
                    start = ((floor((jj-1)/8)-1)*8)+1;
                    flag(ii,start:start+7) = 1; %% prev = 1
                else
                    start = 512-7;
                    flag(ii-1,start:start+7) = 1; %% prev in prev_line = 1
                end
                % flag next = 1
                if floor((jj-1)/8) ~= 63 %% current block ~= last
                    start = ((floor((jj-1)/8)+1)*8)+1; %% next = 1
                    flag(ii,start:start+7) = 1;
                else
                    start = 1;
                    flag(ii+1,start:start+7) = 1; %% next in next_line = 1
                end
            else
                % flag current = 0
                start = (floor((jj-1)/8)*8)+1;
                flag(ii,start:start+7) = 0;
                % flag next = 1
                if (floor((jj-1)/8)*8)+1 ~= 505 %% current block ~= last
                    start = ((floor((jj-1)/8)+1)*8)+1; %% next = 1
                    flag(ii,start:start+7) = 1;
                else
                    start = 1;
                    flag(ii+1,start:start+7) = 1; %% next in next_line = 1
                end
            end
            jj = (floor((jj-1)/8)+1)*8; %% set jj to end of current block
        end
        jj = jj + 1;
    end
end

% correcting [ flag error flag error flag ] to [flag error blank error flag]
for ii = 1:512
    for jj = 1:512
        if mod(jj,8) == 0
           if(jj >= 24)
               if any(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])...
                       && all(flag(ii,jj-15:jj-8) == [1,1,1,1,1,1,1,1])...
                       && any(error(ii,jj-23:jj-16) == [1,1,1,1,1,1,1,1])
                   flag(ii,jj-15:jj-8) = 0;
               end
           elseif(ii~=1 && jj >= 16)
               if any(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])...
                       && all(flag(ii,jj-15:jj-8) == [1,1,1,1,1,1,1,1])...
                       && any(error(ii-1,505:512) == [1,1,1,1,1,1,1,1])
                   flag(ii,jj-15:jj-8) = 0;
               end
           elseif(ii~=1 && jj >= 8)
               if any(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])...
                       && all(flag(ii-1,505:512) == [1,1,1,1,1,1,1,1])...
                       && any(error(ii-1,497:504) == [1,1,1,1,1,1,1,1])
                   flag(ii-1,505:512) = 0;
               end
           end
        end
    end
end

% figure('Name','Flag');
% imshow((flag),[]);

figure('Name','Flag + Error');
imshow((flag+(error*2)),[]);

% figure('Name','Error');
% imshow((error),[]);

%% Data Embedding %%


data = inputdlg({'Enter the Data to be Embedded'},'Embedding Data',[5,35],{'This is the data to be embedded'});
% data = inputdlg({'Enter the Data to be Embedded'},'Embedding Data',[5,35],{data});
data = (cell2mat(data));
data = [data '..........'];

M = [data zeros(1,floor(512*512/8) - numel(data))];
% M = data;

Kw = inputdlg({'Enter the Key for Encrypting Word'},'Word Encryption Key',[1,35],{'123456'});
Kw = str2double((cell2mat(Kw)));
seed = Kw;
rng(seed,'twister');
S = randi(255,[1,numel(M)]);

M = double(M);
Me = bitxor(M,S);
Me = Me';
Me = dec2bin(Me);
Me = Me';
Me = reshape(Me,[1,numel(Me)]);
Me = double(Me);
Me = Me - 48;

check = 0;
empty = ones(512);
for ii = 1:512
    for jj = 1:512
        if flag(ii,jj) == 1
            check = check+1;
        end
        if check > 0
            empty(ii,jj) = 0;
        end
        if check == 16
            check = 0;
        end
    end
end

empty(1,1:16) = 0;
% empty(1,:) = 0;
% empty(:,1) = 0;

% figure('Name','Empty spaces available for hiding','NumberTitle','off');
% imshow(uint8(empty),[]);

Iew = Ie;
m = 1;
for ii = 1:512
    for jj = 1:512
        if ii == 1 && jj <= 8
            continue;
        end
        if empty(ii,jj) == 1
            Iew(ii,jj) = bitset(Iew(ii,jj),8,Me(m));
            m = m + 1;
            if mod(jj,8) == 0
                if all(bitget(Iew(ii,jj-7:jj),8) == [1,1,1,1,1,1,1,1])
                    Iew(ii,jj) = bitset(Iew(ii,jj),8,0);
                end
            end
        else
            Iew(ii,jj) = bitset(Iew(ii,jj),8,(flag(ii,jj) | error(ii,jj)));
        end
    end
end

figure('Name','Encrypted Image with hidden word','NumberTitle','off');
imshow(uint8(Iew));

%% Extraction %%
%% Message Extraction %%

to_check = ones(512);
error_check = zeros(512);
error_extract = zeros(512);
error = zeros(512);
MSBs = bitget(Iew,8);
check = 0;

for ii = 1:512
    for jj = 1:8:512
        if ~all(MSBs(ii,jj:jj+7) == [1,1,1,1,1,1,1,1]) && check == 1
            error(ii,jj:jj+7) = MSBs(ii,jj:jj+7);
        end
        if all(MSBs(ii,jj:jj+7) == [1,1,1,1,1,1,1,1])
            check = check + 1;
        end
        if check > 0
            to_check(ii,jj:jj+7) = 0;
        end
        if check == 2
            check = 0;
        end
    end
end

to_check(1,1:16) = 0;

% figure('Name','Empty spaces being checked','NumberTitle','off');
% imshow(uint8(to_check),[]);
% figure('Name','Difference between Empty spaces available and checked','NumberTitle','off');
% imshow(uint8(empty ~= to_check),[]);

m=0;
Me = zeros(1,512*512);
for ii = 1:512
    for jj = 1:512
        if to_check(ii,jj) == 1
%             Me(m+1) = bitget(Iew(ii,jj),8);
            Me(m+1) = MSBs(ii,jj);
            m = m + 1;
        end
    end
end

Kw = inputdlg({'Enter the Key for Decrypting Word'},'Word Decryption Key',[1,35],{'123456'});
Kw = str2double((cell2mat(Kw)));
seed = Kw;
rng(seed,'twister');
S = randi(255,[1,numel(Me)/8]);

Me = Me + 48;
Me = char(Me);
Me = reshape(Me,[8,numel(Me)/8]);
Me = Me';
Me = bin2dec(Me);
Me = (Me');
Md = bitxor(Me,S);
Md = char(Md);

Md = strsplit(Md,'...');
Md = cell2mat(Md(1));
fprintf(Md);
fprintf('\n\n');
msgbox(Md,'Decoded Message');

%% Image Extraction %%

Ke = inputdlg({'Enter the Key for Decrypting Image'},'Image Decryption Key',[1,35],{'123456789'});
Ke = str2double((cell2mat(Ke)));
seed = Ke;
rng(seed,'twister');
S = randi(255,512);
Id = bitxor(S,Iew);
% figure('Name','Decoded Image','NumberTitle','off');
% imshow(uint8(Id));

% z = 1;
% eqdec(:,1) = [0;0];
for ii = 1:512
    for jj  = 1:512
        Id0m = bitset(Id(ii,jj),8,0);
        Id1m = bitset(Id(ii,jj),8,1);
        if ii == 1 && jj <= 8
            continue;
        elseif ii == 1 && jj ~= 1
            delta0 = abs(Id0m-Id(ii,jj-1));
            delta1 = abs(Id1m-Id(ii,jj-1));
        elseif ii ~= 1 && jj == 1
            delta0 = abs(Id0m-Id(ii-1,jj));
            delta1 = abs(Id1m-Id(ii-1,jj));
        else
            delta0 = abs(Id0m-floor((Id(ii,jj-1)+Id(ii-1,jj))/2));
            delta1 = abs(Id1m-floor((Id(ii,jj-1)+Id(ii-1,jj))/2));
        end
        if delta0 < delta1
            Id(ii,jj) = bitset(Id(ii,jj),8,error(ii,jj));
%             Id(ii,jj) = bitset(Id(ii,jj),8,0);
        elseif delta0 == delta1
            Id(ii,jj) = bitset(Id(ii,jj),8,error(ii,jj));
%             eqdec(:,z) = [ii;jj];
%             z = z+1;
%             fprintf('equal predict [%i %i] \n',ii,jj);
        else
            Id(ii,jj) = bitset(Id(ii,jj),8,~error(ii,jj));
%             Id(ii,jj) = bitset(Id(ii,jj),8,1);
        end
    end
end

figure('Name','Corrected Decoded Image','NumberTitle','off');
imshow(uint8(Id));

figure('Name','Relative Id vs I Uncorrected')
imshow(uint8(abs(Id-I)));

figure('Name','Absolute Id vs I Uncorrected')
imshow((Id~=I),[]);

% [row,col] = find(Id'~=In');
% fprintf('[%i %i]\n', [col,row].');

% [lic,loc] = ismember(eqdec',eqpre','rows');
% loc;


%% PSNR and SSIM

img_processed = Id;
img_original = I;
psnrval = psnr(img_processed,img_original,255);
ssimval = ssim(img_processed,img_original);
fprintf('psnr = %4.2f dB \nssim = %1.5f \n',psnrval,ssimval);
bpp = ((numel(find(empty ==1)))/(512*512));
fprintf('Max bits per pixel = %1.4f \n',bpp);
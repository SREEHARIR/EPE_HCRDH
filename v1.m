clc;
close all;
clear;

Ke = 1234567898;
Kw = 12345;

data = ['This is the message to be embedded. '...
    'This is the continuation of the same message.\n'...
    'This is a really long message to check the amount of problems in the word hiding. \n'...
    'The size of the message is limited to around 512*512/8 characters in the string.\n'...
    'This string should be fully decipherable and should be shown in the command Window.\n\n'];
data = [data data data data data data data data data data data data data data data];

% img = imread('./images/hands1.jpg');
% img = imread('./images/satellite.png');
img = imread('./images/airplane.pgm');

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

pred = zeros(512);
In=zeros(512);
error = zeros(512);

for ii = 1:512
    for jj = 1:512
%         if ii == 1 || jj == 1
%             pred(ii,jj) = I(ii,jj);
%         else
%             pred(ii,jj) = floor((pred(ii-1,jj) + pred(ii,jj-1))/2);
% %             pred(ii,jj) = floor(I(ii-1,jj) + I(ii,jj-1)/2);
%         end
        if (ii == 1) || (jj == 1)
            pred(ii,jj) = I(ii,jj);
            if(jj ~= 1)
                pred(ii,jj) = I(ii,jj-1);
            elseif ii ~= 1
                pred(ii,jj) = I(ii,1,jj);
            end
        else
            if abs(I(ii-1,jj)-I(ii,jj)) < abs(I(ii,jj-1)-I(ii,jj))
                pred(ii,jj) = I(ii-1,jj);
            else
                pred(ii,jj) = I(ii,jj-1);
            end
        end
        delta = abs(pred(ii,jj) - I(ii,jj));
        deltainv = abs(pred(ii,jj) - Iinv(ii,jj));
        if(delta >= deltainv)
            In(ii,jj) = I(ii,jj);
            error(ii,jj) = 1;
        else
            In(ii,jj) = I(ii,jj);
        end
        pred(ii,jj) = In(ii,jj);
%         I(ii,jj) = In(ii,jj);
    end
end

% figure('Name','Prediction Corrected vs Original');
% imshow(In~=I,[]);

% figure('Name','Error');
% imshow((error),[]);
%% Encryption %%

seed = Ke;
rng(seed,'twister');
S = randi(255,512);

Ie = bitxor(S,In);
% figure('Name','Encrypted Image','NumberTitle','off');
% imshow(uint8(Ie));

%% Adapting to avoid prediction error

flag = zeros([512 512]);
for ii = 2:512
    for jj = 1:512
        if error(ii,jj) == 1
            if flag(ii,jj) == 0
%                 flag_previous = 1;
                for k = 1:8
                    if floor((jj-1)/8) ~= 0 %% current_block ~= 0
                        flag(ii,((floor((jj-1)/8)-1)*8)+k) = 1; %% prev = 1
                    else
                        flag(ii-1,512-8+k) = 1; %% prev in prev_line = 1
                    end
                end
%                 flag_next = 1;
                for k = 1:8
                    if floor((jj-1)/8) ~= 63 %% current block ~= last
                        flag(ii,((floor((jj-1)/8)+1)*8)+k) = 1; %% next = 1
                    else
                        flag(ii+1,0+k) = 1; %% next in next_line = 1
                    end
                end
            else
%                 flag_current = 0;
                for k = 1:8
                    flag(ii,(floor((jj-1)/8)*8)+k) = 0; %% current = 0
                end
%                 flag_next = 1;
                for k = 1:8
                    if floor((jj-1)/8) ~= 63 %% current_block ~= last
                        flag(ii,((floor((jj-1)/8)+1)*8)+k) = 1; %% next = 1
                    else
                        flag(ii+1,0+k) = 1; %% next in next_line = 1
                    end
                end
            end
            jj = (floor((jj-1)/8)+1)*8; %% set jj to end of current block
        end
    end
end

% figure('Name','Flag');
% imshow((flag),[]);
% figure('Name','Flag + Error');
% imshow((flag+error),[]);
% figure('Name','Error');
% imshow((error),[]);

%% Data Embedding %%

M = [data zeros(1,floor(512*512/8) - numel(data))];
% M = data;
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

Iew = Ie;
m = 0;
check = 0;
empty = zeros(512);
track = 0;
for ii = 1:512
    for jj = 1:512
        if flag(ii,jj) == 1
            check = check + 1;
        end
        if check == 0 && ii ~= 1 && jj ~= 1
            Iew(ii,jj) = bitset(Iew(ii,jj),8,Me(m+1));
            m = m + 1;
            empty(ii,jj) = 1;
            if mod(jj,8) == 0
                for k = 0:7
                    if bitget(Iew(ii,jj-k),8) == 1
                        track = track + 1;
                    end
                end
                if track == 8
                    m = m - 8;
%                     Iew(ii,jj:512) = bitset(Iew(ii,jj:512),8,1);
                    %set next 8 to zero.
                    %If next 8 is already zero, set to 1
                    if jj == 512
                        ii = ii + 1;
                        jj = 1;
                        Iew(ii,jj:jj+7) = bitset(Iew(ii,jj:jj+7),8,~bitget(Iew(ii,jj:jj+7),8));
                        jj = 8;
                    else
                        jj = jj + 1;
                        Iew(ii,jj:jj+7) = bitset(Iew(ii,jj:jj+7),8,~bitget(Iew(ii,jj:jj+7),8));
                        jj = jj + 7;
                    end
                end
                track = 0;
            end
        elseif check ~= 0 && ii ~=1 && jj ~= 1
            Iew(ii,jj) = bitset(Iew(ii,jj),8,(flag(ii,jj)+error(ii,jj)));
        end
        if check == 16
            check = 0;
        end
    end
end

figure('Name','Empty spaces available for hiding','NumberTitle','off');
imshow(uint8(empty),[]);

% figure('Name','Encrypted Image with hidden word','NumberTitle','off');
% imshow(uint8(Iew));

%% Extraction %%
%% Message Extraction %%

m = 0;
check = 0;
Me = zeros(1,512*512);
checking = zeros(512);
for ii = 2:512
%     fprintf('%d \n',ii);
    for jj = 1:512
        if ii ~= 1 && jj ~= 1 && check == 0
            Me(m+1) = bitget(Iew(ii,jj),8);
            m = m + 1;
            checking(ii,jj) = 1;
        end
        if mod(jj,8) == 0
            for k = 0:7
                if bitget(Iew(ii,jj-k),8) == 1
                    check = check + 1;
                end
            end
            if check == 8
                Me(m-7:m) = zeros([1,8]);
                m = m - 8;
            elseif check > 0 && check < 8
                check = 0;
            elseif check == 16
                check = 0;
            elseif check > 8 && check < 16
                check = 8;
            end
        end
    end
end

figure('Name','Spaces being checked','NumberTitle','off');
imshow(uint8(checking),[]);

Me = Me(1:261120); %% floor((512*512-512-511-numel(find(flag==1))-numel(find(error==1)))/8)*8

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


fprintf(Md);
fprintf('\n');

%% Image Extraction %%

seed = Ke;
rng(seed,'twister');
S = randi(255,512);
Id = bitxor(S,Iew);
% figure('Name','Decoded Image','NumberTitle','off');
% imshow(uint8(Id));

% figure;
for ii = 1:512
    for jj  = 1:512
%         predictor = floor((Id(ii-1,jj) + Id(ii,jj-1))/2);
        
        if (ii == 1) || (jj == 1)
            predictor = Id(ii,jj);
            if jj ~= 1
                predictor = Id(ii,jj-1);
            elseif ii ~= 1
                predictor = Id(ii-1,jj);
            end
        else
            if abs(Id(ii-1,jj)-Id(ii,jj)) < abs(Id(ii,jj-1)-Id(ii,jj))
                predictor = Id(ii-1,jj);
            else
                predictor = Id(ii,jj-1);
            end
        end
        
        Id0m = bitset(Id(ii,jj),8,0);
        Id1m = bitset(Id(ii,jj),8,1);
        delta0 = abs(predictor - Id0m);
        delta1 = abs(predictor - Id1m);
        if delta0 < delta1
            Id(ii,jj) = bitset(Id(ii,jj),8,0);
%             Id(ii,jj) = bitset(Id(ii,jj),8,error(ii,jj));
        else
            Id(ii,jj) = bitset(Id(ii,jj),8,1);
%             Id(ii,jj) = bitset(Id(ii,jj),8,~error(ii,jj));
        end
    end
%     imshow(uint8(Id));
end


% figure('Name','Predicted');
% imshow(uint8(pred));

figure('Name','Corrected Decoded Image','NumberTitle','off');
imshow(uint8(Id));

figure('Name','Id vs I Uncorrected')
imshow(((Id~=I)),[]);


%% PSNR and SSIM

img_processed = Id;
img_original = I;
psnrval = psnr(img_processed,img_original,255);
ssimval = ssim(img_processed,img_original);
fprintf('psnr = %4.2f dB \nssim = %1.5f \n',psnrval,ssimval);

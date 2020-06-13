clc;
close all;
clear;

fprintf('program Start');

Ke = 12345678;
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

% pred = zeros(512);
In=zeros(512);
error = zeros(512);

for ii = 2:512
    for jj = 2:512
%         if (ii == 1) || (jj == 1)
%             predictor = I(ii,jj);
% %             if(jj ~= 1)
% %                 pred(ii,jj) = I(ii,jj-1);
% %             elseif ii ~= 1
% %                 pred(ii,jj) = I(ii,1,jj);
% %             end
%         else
%             predictor = floor((I(ii-1,jj) + I(ii,jj-1))/2);
%             
% %             deltau = abs(I(ii-1,jj)-I(ii,jj));
% %             deltal = abs(I(ii,jj-1)-I(ii,jj));
% % %             if deltal == deltau
% % %                 I(ii,jj) = I(ii,jj) - 1;
% % %             end
% %             if deltau < deltal
% %                 predictor = I(ii-1,jj);
% %             else
% %                 predictor = I(ii,jj-1);
% %             end
%         end
%         delta = abs(predictor - I(ii,jj));
%         deltainv = abs(predictor - Iinv(ii,jj));
        delta = min(abs(I(ii,jj)-I(ii,jj-1)),abs(I(ii,jj)-I(ii-1,jj)));
        deltainv = min(abs(Iinv(ii,jj)-I(ii,jj-1)),abs(Iinv(ii,jj)-I(ii-1,jj)));
        if(delta >= deltainv)
            In(ii,jj) = I(ii,jj);
            error(ii,jj) = 1;
%         elseif delta == deltainv
%             In(ii,jj) = I(ii,jj)+((predictor - I(ii,jj))/delta);
        else
            In(ii,jj) = I(ii,jj);
        end
%         pred(ii,jj) = In(ii,jj);
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

%% Create Flag bits to surround error bits

flag = zeros([512 512]);
% figure;
for ii = 1:512
    jj = 1;
    while jj <= 512
%         imshow(flag+error,[]);
% display(jj);
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


for ii = 1:512
    for jj = 1:512
        if mod(jj,8) == 0
           if(jj >= 24)
               if any(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])...
                       && all(flag(ii,jj-15:jj-8) == [1,1,1,1,1,1,1,1])...
                       && any(error(ii,jj-23:jj-16) == [1,1,1,1,1,1,1,1])
                   flag(ii,jj-15:jj-8) = 0;
               end
           elseif(ii~=0 && jj >= 16)
               if any(error(ii,jj-7:jj) == [1,1,1,1,1,1,1,1])...
                       && all(flag(ii,jj-15:jj-8) == [1,1,1,1,1,1,1,1])...
                       && any(error(ii-1,505:512) == [1,1,1,1,1,1,1,1])
                   flag(ii,jj-15:jj-8) = 0;
               end
           elseif(ii~=0 && jj >= 8)
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

% figure('Name','Flag + Error');
% imshow((flag+(error*2)),[]);

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

empty(1,:) = 0;
empty(:,1) = 0;

% figure('Name','Empty spaces available for hiding','NumberTitle','off');
% imshow(uint8(empty),[]);

for ii = 1:512
    jj = 1;
    while jj <= 512
        if empty(ii,jj) == 1
            Iew(ii,jj) = bitset(Iew(ii,jj),8,Me(m+1));
            m = m + 1;
%             fprintf('m = %i\n',m);
            if mod(jj,8) == 0
%                 imshow(uint8(bitget(Iew,8)),[]);
                if all(bitget(Iew(ii,jj-7:jj),8) == [1,1,1,1,1,1,1,1])
%                     fprintf('m = %i\n',m);
%                     if jj ~= 512
%                         if all(flag(ii,jj+1:jj+8) ~= [1,1,1,1,1,1,1,1]) %% wont work because jj+1 to jj+8 has not yet been assigned ??? No... %% will work but should also edit flag to remove jj+1 to jj+8, not just Iew No... maybe this is all correct
%                             Iew(ii,jj+1:jj+8) = bitset(Iew(ii,jj+1:jj+8),8,1);
%                         elseif ii~= 512
%                             Iew(ii,jj+1:jj+8) = bitset(Iew(ii,jj+1:jj+8),8,0);
%                             fprintf('removing flag in 1st else');
%                         end
%                     else
%                         if ii~=512 && all(flag(ii+1,1:8) ~= [1,1,1,1,1,1,1,1])
%                             Iew(ii+1,1:8) = bitset(Iew(ii+1,1:8),8,1);
%                         elseif ii~=512
%                             Iew(ii+1,1:8) = bitset(Iew(ii+1,1:8),8,0);
%                             fprintf('removing flag in 2nd else');
%                         end
%                     end
%                     jj = jj + 8;
%                     m = m - 8;
%                     m = m - 1;
                    Iew(ii,jj) = bitset(Iew(ii,jj),8,0);
                end
            end
        else
            Iew(ii,jj) = bitset(Iew(ii,jj),8,(flag(ii,jj) | error(ii,jj)));
        end
        jj = jj + 1;
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
%             error_extract(ii,jj:jj+7) = MSBs(ii,jj:jj+7);
            error(ii,jj:jj+7) = MSBs(ii,jj:jj+7);
        end
        if all(MSBs(ii,jj:jj+7) == [1,1,1,1,1,1,1,1])
            check = check + 1;
        end
        if check > 0
            to_check(ii,jj:jj+7) = 0;
        end
%         if check == 1
%             error_check(ii,jj:jj+7) = 1;
%         end
        if check == 2
            check = 0;
        end
    end
end

m=0;
Me = zeros(1,512*512);
for ii = 2:512
    for jj = 2:512
        if to_check(ii,jj) == 1
            Me(m+1) = bitget(Iew(ii,jj),8);
            m = m + 1;
        end
    end
end

% figure('Name','Spaces being checked for message','NumberTitle','off');
% imshow(uint8(to_check),[]);
% 
% figure('Name','Spaces checked for error stored','NumberTitle','off');
% imshow(uint8(error_check),[]);
% 
% figure('Name','Extracted Error','NumberTitle','off');
% imshow(uint8(error_extract),[]);

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
% figure('Name','Id vs I Uncorrected')
for ii = 2:512
%     imshow(((Id~=I)),[]);
    for jj  = 2:512
%         predictor = floor((Id(ii-1,jj) + Id(ii,jj-1))/2);
        
%         if (ii == 1) || (jj == 1)
%             predictor = Id(ii,jj);
% %             if jj ~= 1
% %                 predictor = Id(ii,jj-1);
% %             elseif ii ~= 1
% %                 predictor = Id(ii-1,jj);
% %             end
%         else
%             deltau = abs(Id(ii-1,jj)-Id(ii,jj));
%             deltal = abs(Id(ii,jj-1)-Id(ii,jj));
%             
%             deltau_c = abs(I(ii-1,jj)-I(ii,jj));
%             deltal_c = abs(I(ii,jj-1)-I(ii,jj));
%             
% %             if deltau ~= deltau_c
% %                 fprintf('\n ii-%i,jj-%i    u_d-%i u_o-%i    l_d-%i l_o-%i    I_d-%i I_o-%i',ii,jj,Id(ii-1,jj),I(ii-1,jj),Id(ii,jj-1),I(ii,jj-1),Id(ii,jj),I(ii,jj));
% %                 fprintf('\n');
% %             end
% %             if deltal ~= deltal_c
% %                 fprintf('\n ii-%i,jj-%i    u_d-%i u_o-%i    l_d-%i l_o-%i    I_d-%i I_o-%i',ii,jj,Id(ii-1,jj),I(ii-1,jj),Id(ii,jj-1),I(ii,jj-1),Id(ii,jj),I(ii,jj));
% %                 fprintf('\n');
% %             end
%             
%             if deltau < deltal
%                 predictor = Id(ii-1,jj);
%             else
%                 predictor = Id(ii,jj-1);
%             end
%         end
        
        Id0m = bitset(Id(ii,jj),8,0);
        Id1m = bitset(Id(ii,jj),8,1);
%         delta0 = abs(predictor - Id0m);
%         delta1 = abs(predictor - Id1m);
        delta0 = min(abs(Id0m-I(ii,jj-1)),abs(Id0m-I(ii-1,jj)));
        delta1 = min(abs(Id1m-I(ii,jj-1)),abs(Id1m-I(ii-1,jj)));
        if delta0 < delta1
            Id(ii,jj) = bitset(Id(ii,jj),8,error(ii,jj));
%             Id(ii,jj) = Id0m;
        else
            Id(ii,jj) = bitset(Id(ii,jj),8,~error(ii,jj));
%             Id(ii,jj) = Id1m;
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

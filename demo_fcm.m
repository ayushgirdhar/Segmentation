clear;close all;
AC=0;
format short

im=imread('6.png');
im=imresize(im,[128,128]);
he=im;
he1=he;
im12=he1./255;
figure(1)
imshow(he), title('H&E image');
he1= imnoise(he,'salt & pepper', 0.05);
cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

[center,member,klm]=fcm(ab,3);
[center,cidx]=sort(center);
member=member';
member=member(:,cidx);
[maxmember,label]=max(member,[],2);
fismat = genfis2(ab,label,[0.5 0.25 0.3]);
pixel_labels = reshape(label,nrows,ncols);

% figure(2)
% imshow(pixel_labels,[]), title('image labeled by cluster index');
segmented_images = cell(1,5);
rgb_label = repmat(pixel_labels,[1 1 3]);
%%


for k=1:3
    for i=1:nrows
        for j=1:ncols
            if(rgb_label(i,j,k)==1)
                im12(i,j,k)=im12(i,j,k)+1;
            elseif(rgb_label(i,j,k)==2)
                im12(i,j,k)=im12(i,j,k)+2;
                
            elseif(rgb_label(i,j,k)==3)
                im12(i,j,k)=im12(i,j,k)+3;
                
            elseif(rgb_label(i,j,k)==4)
                im12(i,j,k)=im12(i,j,k)+4;
                
            elseif(rgb_label(i,j,k)==5)
                im12(i,j,k)=im12(i,j,k)+5;
            
            elseif(rgb_label(i,j,k)==6)
                im12(i,j,k)=im12(i,j,k)+6;
           
            elseif(rgb_label(i,j,k)==7)
                im12(i,j,k)=im12(i,j,k)+7;
            
            elseif(rgb_label(i,j,k)==8)
                im12(i,j,k)=im12(i,j,k)+8;
            
           
            elseif(rgb_label(i,j,k)==9)
                im12(i,j,k)=im12(i,j,k)+9;
        
             
            elseif(rgb_label(i,j,k)==10)
                im12(i,j,k)=im12(i,j,k)+10;
            end
            
        end
    end
end


for k = 1:10
    color = he;
    color(rgb_label ~= k) = 0;
    if k==1
    ind_1=find(color(rgb_label == k));
    ind_11=find(color == 0);
    elseif k==2
    ind_2=find(color(rgb_label == k));
    ind_22=find(color == 0);
    elseif k==3
    ind_3=find(color(rgb_label == k));
    ind_33=find(color== 0);
    elseif k==4
    ind_4=find(color(rgb_label == k));
    ind_44=find(color == 0);
    elseif k==5
    ind_5=find(color(rgb_label == k));
    ind_55=find(color == 0);
    elseif k==6
    ind_2=find(color(rgb_label == k));
    ind_22=find(color == 0);
    elseif k==7
    ind_3=find(color(rgb_label == k));
    ind_33=find(color== 0);
    elseif k==8
    ind_4=find(color(rgb_label == k));
    ind_44=find(color == 0);
    elseif k==9
    ind_5=find(color(rgb_label == k));
    ind_55=find(color == 0);
    elseif k==10
    ind_5=find(color(rgb_label == k));
    ind_55=find(color == 0);
    
    end
    
    segmented_images{k} = color;
end

%% new addition
new_im12=im12;new_im13=im12;new_im14=im12;new_im15=im12;new_im16=im12;

new_im12(ind_11)=0;new_im13(ind_22)=0;new_im14(ind_33)=0;new_im15(ind_44)=0;new_im16(ind_55)=0;

sw=size(he);
total_pix=sw(1)*sw(2)*sw(3);
format long
total_p_1=(numel(ind_1)/total_pix)*100;
total_p_2=(numel(ind_2)/total_pix)*100;
total_p_3=(numel(ind_3)/total_pix)*100;
total_p_4=(numel(ind_4)/total_pix)*100;
total_p_5=(numel(ind_5)/total_pix)*100;
figure(2);
colormap(gray);
imagesc(im12(:,:,1));

figure(3);
subplot(3,2,1)
 imshow(he), title('H&E image');
 subplot(3,2,2)
 str1=sprintf('objects in cluster (%f)',total_p_1);
imshow(segmented_images{1}), title(str1);
subplot(3,2,3)
str1=sprintf('objects in cluster (%f)',total_p_2);
imshow(segmented_images{2}), title(str1);
subplot(3,2,4)
str1=sprintf('objects in cluster (%f)',total_p_3);
imshow(segmented_images{3}), title(str1);
subplot(3,2,5)
str1=sprintf('objects in cluster (%f)',total_p_4);
imshow(segmented_images{4}), title(str1);
subplot(3,2,6)
str1=sprintf('objects in cluster (%f)',total_p_5);
imshow(segmented_images{5}), title(str1);

% 
% mean_cluster_value = mean(center,2);
% [tmp, idx] = sort(mean_cluster_value);
% pink_cluster_num = idx(4);
%  nrols=23000;
% L = lab_he(:,:,1);
% pink_idx = find(pixel_labels == pink_cluster_num);
% L_pink = L(pink_idx);
% is_light_pink = im2bw(L_pink,graythresh(L_pink));
% wave_labels = repmat(uint8(0),[nrows ncols]);
% wave_labels(pink_idx(is_light_pink==false)) = 1;
% wave_labels = repmat(wave_labels,[1 1 3]);
% pink_wave = he;
% pink_wave(wave_labels ~= 1) = 0;
% figure(4);
% imshow(pink_wave), title('pink wave');

figure();
imshow(he1);title('noisy image');

lab_he1 = applycform(he1,cform);
ab1 = double(lab_he1(:,:,2:3));
nrows = size(ab1,1);
ncols = size(ab1,2);
ab1 = reshape(ab1,nrows*ncols,2);
 
nColors = 5;
% repeat the clustering 3 times to avoid local minima
[center1,member1]=fcm(ab1,5);
[center1,cidx1]=sort(center1);
member1=member1';
member1=member1(:,cidx1);
[maxmember1,label1]=max(member1,[],2);

label1 = bestMap(label,label1);
AC = length(find(label == label1))/length(label);
fprintf('Accuracy :: %f\n\n',AC);
pixel_labels1 = reshape(label1,nrows,ncols);

figure();
subplot(211);
imshow(pixel_labels,[]), title('image labeled by cluster index without noise');
subplot(212);
imshow(pixel_labels1,[]), title('image labeled by cluster index in noisy');

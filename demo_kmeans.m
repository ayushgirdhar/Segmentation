close all; clear all;clc;

format short


he = imread('7.png');
he1=double(he);
im12=he1./255;
he1= imnoise(he,'salt & pepper', 0.05);
cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 5;

cluster_center_init=[130.201,131.388;115.134,159.846;137.952,157.969;172.768,132.083;150.493,177.837];
clus_init(:,:,1)=cluster_center_init;clus_init(:,:,6)=cluster_center_init;
clus_init(:,:,2)=cluster_center_init;clus_init(:,:,7)=cluster_center_init;
clus_init(:,:,3)=cluster_center_init;clus_init(:,:,8)=cluster_center_init;
clus_init(:,:,4)=cluster_center_init;clus_init(:,:,9)=cluster_center_init;
clus_init(:,:,5)=cluster_center_init;clus_init(:,:,10)=cluster_center_init;

% repeat the clustering 3 times to avoid local minima
[cluster_idx ,cluster_center] = kmeans(ab,5,'distance','sqEuclidean', ...
                                      'replicate',10,'start',clus_init,'Display','Iter');

pixel_labels = reshape(cluster_idx,nrows,ncols);
figure(2)
imshow(pixel_labels,[]), title('image labeled by cluster index');
rgb_label = repmat(pixel_labels,[1 1 3]);
%% new addition
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
%%
segmented_images = cell(1,5);

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
sw=size(he);
total_pix=sw(1)*sw(2)*sw(3);
format long
total_p_1=(numel(ind_1)/total_pix)*100;
total_p_2=(numel(ind_2)/total_pix)*100;
total_p_3=(numel(ind_3)/total_pix)*100;
total_p_4=(numel(ind_4)/total_pix)*100;
total_p_5=(numel(ind_5)/total_pix)*100;

%% new addition
new_im12=im12;new_im13=im12;new_im14=im12;new_im15=im12;new_im16=im12;

new_im12(ind_11)=0;new_im13(ind_22)=0;new_im14(ind_33)=0;new_im15(ind_44)=0;new_im16(ind_55)=0;
%%


% figure(4);
% subplot(3,2,1)
%  imshow(he), title('H&E image');
%  subplot(3,2,2)
%  str1=sprintf('objects in cluster');
% imshow(segmented_images{6}), title(str1);
% subplot(3,2,3)
% str1=sprintf('objects in cluster');
% imshow(segmented_images{7}), title(str1);
% subplot(3,2,4)
% str1=sprintf('objects in cluster');
% imshow(segmented_images{8}), title(str1);
% subplot(3,2,5)
% str1=sprintf('objects in cluster');
% imshow(segmented_images{9}), title(str1);
% subplot(3,2,6)
% str1=sprintf('objects in cluster');
% imshow(segmented_images{10}), title(str1);


lab_he = applycform(he1,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab1 = reshape(ab,nrows*ncols,2);
 
nColors = 5;
cluster_center_init=[130.201,131.388;115.134,159.846;137.952,157.969;172.768,132.083;150.493,177.837];
clus_init(:,:,1)=cluster_center_init;clus_init(:,:,6)=cluster_center_init;
clus_init(:,:,2)=cluster_center_init;clus_init(:,:,7)=cluster_center_init;
clus_init(:,:,3)=cluster_center_init;clus_init(:,:,8)=cluster_center_init;
clus_init(:,:,4)=cluster_center_init;clus_init(:,:,9)=cluster_center_init;
clus_init(:,:,5)=cluster_center_init;clus_init(:,:,10)=cluster_center_init;

% repeat the clustering 3 times to avoid local minima
[cluster_idx1 ,cluster_center1] = kmeans(ab1,5,'distance','sqEuclidean', ...
                                      'replicate',10,'start',clus_init,'Display','Iter');


                                  
cluster_idx1 = bestMap(cluster_idx,cluster_idx1);
                                  
% ct1=zeros(size(cluster_idx1,1),size(cluster_idx1,2));
% ct1(find(cluster_idx1==3))=1;
% ct1(find(cluster_idx1==1))=3;
% ct1(find(cluster_idx1==2))=2;
                                  % 
% mean_cluster_value = mean(cluster_center,2);
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

AC = length(find(cluster_idx == cluster_idx1))/length(cluster_idx);


figure(2);
colormap(gray);
imagesc(im12(:,:,1));

figure(3);
subplot(3,2,1)
 imshow(he), title('H&E image');
 subplot(3,2,2)
 str1=sprintf('objects in cluster');
imshow(segmented_images{1}), title(str1);
subplot(3,2,3)
str1=sprintf('objects in cluster');
imshow(segmented_images{2}), title(str1);
subplot(3,2,4)
str1=sprintf('objects in cluster');
imshow(segmented_images{3}), title(str1);
subplot(3,2,5)
str1=sprintf('objects in cluster');
imshow(segmented_images{4}), title(str1);
subplot(3,2,6)
str1=sprintf('objects in cluster');
imshow(segmented_images{5}), title(str1);

figure();
imshow(he1);title('noisy image');

figure(1);
imshow(he), title('H&E image');



fprintf('Accuracy :: %f\n\n',AC);
pixel_labels1 = reshape(cluster_idx1,nrows,ncols);

figure();
imshow(pixel_labels1,[]), title('image labeled by cluster index in noisy');

figure();
subplot(211);
imshow(pixel_labels,[]), title('image labeled by cluster index without noise');
subplot(212);
imshow(pixel_labels1,[]), title('image labeled by cluster index in noisy');

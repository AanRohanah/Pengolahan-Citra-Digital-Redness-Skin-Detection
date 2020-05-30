citra = imread('merah.jpg');
% ---Kode Program Resize Foto---

[bar, kol] = size(citra);

if (bar > kol)
    maxLength = bar;
    if(maxLength >= 500)
        citra = imresize(citra, [480 NaN]);
    end
else
    maxLength = kol;
    if(maxLength >= 500)
        citra = imresize(citra, [NaN 480]);
    end
end

%--- Kode Program Segmentasi Kulit---
%Segemntasi dengan HSV

[bar, kol, dlm] = size(citra);
citrahsv = rgb2hsv(citra);

hue = citrahsv(:,:,1);
sat = citrahsv(:,:,2);
val = citrahsv(:,:,3);

filterByHS= uint8(zeros(bar, kol, dlm));
for i = 1 : bar
    for j = 1 : kol
        if (hue(i, j) <= 0.25 && sat(i, j) >= 0.15 && sat(i, j) <= 0.9)
            filterByHS(i, j, :) = citra(i, j, :);
        end
    end
end

%--- Kode Program Clustering K-Means ---
%Clustering dengan HSV
hs = double(citrahsv(:,:,1:2));
nbar = size (hs,1);
nkol = size (hs,2);
hs = reshape(hs,nbar*nkol,2);

%Membagi kedalam beberapa cluster
nColors = 3; %banyak clustering
[cluster_idx,cluster_center] = kmeans(hs, nColors,'distance','sqEuclidean','Replicates',3);
pixel_labels = reshape(cluster_idx,nbar,nkol);

%Menampilkan hasil segmentasi
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);
for k = 1:nColors
     color = citra;
     color(rgb_label ~= k ) = 0;
     segmented_images{k} = color;
end
 
%Menghitung cluster terbanyak pada daerah wajah
clusterCount = zeros(nColors);
for i = 1:nbar
     for j = 1:nkol
         if filterByHS(i, j, :) ~= 0
             clusterCount (pixel_labels(i,j)) = clusterCount(pixel_labels(i,j)) + 1;
         end
     end
end
[maxVal, maxClusterIndex] = max(clusterCount);
figure,imshow(segmented_images{maxClusterIndex(1)}),
title(strcat(['Objects in Cluster ',num2str(maxClusterIndex(1))]));

%--- Kode Program Imfill ---
%imfill untuk meminimalisir lubang-lubang pada gambar
citraRgb = rgb2gray(segmented_images{maxClusterIndex(1)});
level = graythresh(citraRgb);
citraImfil = imbinarize(citraRgb);
citraImfil = imfill(citraImfil, 'holes');

%--- Kode Program Mengubah Citra Menjadi Berwarna ---
%--- Ekstraksi Ciri (Metode Redness) ---
%mengubah citra menjadi gambar warna
citraImfilWarna = zeros(bar,kol,3);

for i = 1 : bar
    for j = 1 : kol
        if (citraImfil(i,j) > 0)
            citraImfilWarna(i,j,:)=citra(i,j,:);
        end
    end
end

citraImfilWarna=uint8(citraImfilWarna);

figure,
subplot (1,2,1), imshow(citra),
subplot (1,2,2), imshow(citraImfilWarna);

%--- Kode Program Metode Redness ---
%metode redness
%mengambil nilai RGB per layer
I = im2double(citraImfilWarna);
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);

%mencari redness per pixel
redness = zeros(bar,kol);
for i = 1 : bar
    for j = 1 : kol
        redness(i,j) = max(0,((2*R(i,j))-(G(i,j)+B(i,j)))/R(i,j))^2;
    end
end

%--- Kode Program Penandaan Kemerahan dengan Threshold
%seleksi bagian wajah yang lebih dari threshold dikategorikan kemerahan
threshold = median(redness);
citraMerahWarna = citraImfilWarna;
for i = 1 : bar
    for j = 1 : kol
        if redness(i,j) > threshold
            citraMerahWarna(i,j,:)=[255, 0, 0,];
        end
    end
end
figure,
subplot(1,3,1),imshow(citraImfilWarna), title('Citra Asli');
subplot(1,3,2),imshow(redness), title('Deteksi Redness');
subplot(1,3,3),imshow(citraMerahWarna), title('Menandai Bagian Redness');

figure, imshow(citraImfil), title('Citra Biner');

%--- Kode Program Perbaikan Citra(Filtering) ---
gray = rgb2gray(citraMerahWarna);
gaussian = imgaussfilt(gray);
figure, imshow(gaussian), title('Filter Gaussian');

%--- Kode Program Eliminasi Indeks ---
%--- Ekstraksi ciri (Luas dan Warna) ---
%Eliminasi Indeks
rednessBiner = zeros(bar,kol);
elIndex = zeros(bar, kol);
for i = 1 : bar
    for j = 1 : kol
        if gaussian(i,j) == 76
            rednessBiner(i,j) = citraMerahWarna(i,j);
        if redness(i,j) < 1
            elIndex(i,j) = rednessBiner(i,j);
        end
        end
    end
end

figure,
subplot(1,2,1), imshow(rednessBiner), title('Citra Redness Biner');
subplot(1,2,2), imshow(elIndex), title('Eliminasi Indeks');

%--- Kode Program Ekstraksi Ciri Luas ---
%seleksi luas 
candidate = logical(elIndex);
[labeledCandidate, numberOfCandidates] = bwlabel(candidate, 8);

stats = regionprops(labeledCandidate, 'Area');
allArea = [stats.Area];

meanArea = mean(allArea);
stdArea = std(allArea);

indexBlob = find(allArea >= 91);
ambilBlob = ismember(labeledCandidate, indexBlob);

blobBW = ambilBlob > 0;
[labeledBlob, numberOfBlobs] = bwlabel(blobBW);

numberOfBlobs;

figure, imshow(labeledBlob), title('Eliminasi Luas Area');

%--- Kode Program Mean Intensity RGB ---
%Eliminasi berdasarkan warna RGB
red = citraImfilWarna(:, :, 1);
green = citraImfilWarna(:, :, 2);
blue = citraImfilWarna(:, :, 3);

r = regionprops(labeledBlob, red, 'MeanIntensity');
g = regionprops(labeledBlob, green, 'MeanIntensity');
b = regionprops(labeledBlob, blue, 'MeanIntensity');

fiturR = [r.MeanIntensity]';
fiturG = [g.MeanIntensity]';
fiturB = [b.MeanIntensity]';
fiturRGB = [fiturR fiturG fiturB];

meanR = mean(fiturR);
meanG = mean(fiturG);
meanB = mean(fiturB);

stdR = std(fiturR);
stdG = std(fiturG);
stdB = std(fiturB);

indexKemerahan = [];
for i = 1 : numberOfBlobs
    if(fiturR(i) >= (meanR-stdR*1.32) && fiturR(i) <= (meanR+stdR*1.1) && fiturG(i) >= (meanG-stdG*1.32) && fiturG(i) <= (meanG+stdG*1.32) && fiturB(i) >= (meanB-stdB*1.32) && fiturB(i) <= (meanB+stdB*1.32))
    indexKemerahan = [indexKemerahan i];
    end
end

kemerahanBW = ismember(labeledBlob, indexKemerahan);
figure,
subplot(1, 2, 1), imshow(labeledBlob), title('Eliminasi Candidate')
subplot(1, 2, 2), imshow(kemerahanBW), title('Eliminasi Berdasarkan Warna');

%mengubah kemerahanBW menjadi gambar warna
kemerahanWarna = zeros(bar,kol,3);
for i = 1 : bar
    for j = 1 : kol
        if(kemerahanBW(i,j) > 0)
            kemerahanWarna(i,j,:)=citraImfilWarna(i,j,:);
        end
    end
end

kemerahanWarna=uint8(kemerahanWarna);
figure,
subplot(1, 2, 1), imshow(citra), title('Citra Asli');
subplot(1, 2, 2), imshow(kemerahanWarna), title('kemerahanBW menjadi Warna');

%--- Kode Program Mean Intensity HSV ---
%meanintensity dengan HSV
kemerahanHSV= rgb2hsv(kemerahanWarna);
[kemerahanLabeledBlob, numberKemerahanOfBlobs] = bwlabel(kemerahanBW);

hueKemerahan = kemerahanHSV(:, :, 1);
satKemerahan = kemerahanHSV(:, :, 2);
valKemerahan = kemerahanHSV(:, :, 3);

h = regionprops(kemerahanLabeledBlob, hueKemerahan, 'MeanIntensity');
s = regionprops(kemerahanLabeledBlob, satKemerahan, 'MeanIntensity');
v = regionprops(kemerahanLabeledBlob, valKemerahan, 'MeanIntensity');

fiturH = [h.MeanIntensity]';
fiturS = [s.MeanIntensity]';
fiturV = [v.MeanIntensity]';
fitur = [fiturH fiturS fiturV];

meanH = mean(fiturH);
meanS = mean(fiturS);
meanV = mean(fiturV);

stdH = std(fiturH);
stdS = std(fiturS);
stdV = std(fiturV);

indexKemerahanHSV = [];
for i = 1 : numberKemerahanOfBlobs
    if(fiturH(i) >= (meanH-stdH*1.16) && fiturH(i) <= (meanH+stdH*0.5))
        indexKemerahanHSV = [indexKemerahanHSV i];
    end
end

kemerahanHSV = ismember(kemerahanLabeledBlob, indexKemerahanHSV);

figure,
subplot(1, 2, 1), imshow(kemerahanLabeledBlob), title('Candidate Kemerahan');
subplot(1, 2, 2), imshow(kemerahanHSV), title('Eliminasi Berdasarkan Warna HSV');

%--- Kode Program Marking ---
%penandaan pada citra asli
kemerahanEdge = edge(kemerahanHSV, 'canny');
hasil = citra;

for i = 1 : bar
    for j = 1 : kol
        if kemerahanEdge(i, j) == 1
            hasil(i, j, 1) = 255;
            hasil(i, j, 2) = 0;
            hasil(i, j, 3) = 0;
        end
    end
end

hasil = uint8(hasil);
figure,
subplot(1,2,1),imshow(citra),title('Citra Asli')
subplot(1,2,2),imshow(hasil),title('Penandaan Kemerahan');
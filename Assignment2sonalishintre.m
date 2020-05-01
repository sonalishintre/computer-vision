%===================================================
% Computer Vision Programming Assignment 2
% @sonali shintre
% ID: 7602
%===================================================

% ---------------- Question 1 ------------------------
%inputimage = 'bird.jpg';
inputimage = 'IDPicture.png';

I1 = imread(inputimage);
No1 = figure;
subplot(1,2,1)
imshow(I1)
subplot(1,2,2)
imhist(I1)

% to create intensity image
I2 = I1(:,:,1) * 0.299 + I1(:,:,2) * 0.587 + I1(:,:,3) * 0.114;
No2 = figure;
subplot(1,2,1)
imshow(I2)
subplot(1,2,2)
imhist(I2)

I2_imadjust = imadjust(I2);
No3 = figure;
subplot(1,2,1)
imshow(I2_imadjust)
subplot(1,2,2)
imhist(I2_imadjust)

I2_ihisteq = histeq(I2);
No4 = figure;
subplot(1,2,1)
imshow(I2_ihisteq)
subplot(1,2,2)
imhist(I2_ihisteq)


level = multithresh(I2);
I2_thresh = img_thresh(I2, level);
No5 = figure;
subplot(1,2,1)
imshow(I2_thresh)
subplot(1,2,2)
histogram(I2_thresh)

%=========================================================
% ---------------- question 2 ------------------------
% 1X2 Operator
No6 = figure;

Gx=[-1.0 1.0];
Gy=[-1.0;1.0];

subplot(2,2,1);
imshow(I2);
title('Original image');

gradx=conv2(I2,Gx,'same');
gradx1=abs(gradx);
gradx1=norm(gradx1);
subplot(2,2,2);
imshow(gradx1,[]);
title('1X2 Horizontal');
 
grady=conv2(I2,Gy,'same');
grady1=abs(grady);
grady1=norm(grady1);
subplot(2,2,3);
imshow(grady1,[]);
title('1X2 Vertical');
 
gradxy1=gradx1+grady1;
subplot(2,2,4);
imshow(gradxy1,[]);
title('1X2 Combination');

% Sobel Operator
No7 = figure;

Gx=[-1.0 -2.0 -1.0; 0.0 0.0 0.0; 1.0 2.0 1.0];
Gy=[-1.0 0.0 1.0; -2.0 0.0 2.0; -1.0 0.0 1.0];

subplot(2,2,1);
imshow(I2);
title('Original image');

gradx=conv2(I2,Gx,'same');
gradx2=abs(gradx);
gradx2=norm(gradx2);
subplot(2,2,2);
imshow(gradx2,[]);
title('Sobel Horizontal');
 
grady=conv2(I2,Gy,'same');
grady2=abs(grady);
grady2=norm(grady2);
subplot(2,2,3);
imshow(grady2,[]);
title('Sobel Vertical');
 
gradxy2=gradx2+grady2;
subplot(2,2,4);
imshow(gradxy2,[]);
title('Sobel Combination');

% Sobel substract 1X2
No8 = figure;

subplot(2,2,1);
imshow(I2);
title('Original image');

gradx3=abs(gradx2-gradx1);
subplot(2,2,2);
imshow(gradx3,[]);
title('Substract Horizontal');
 
grady3=abs(grady2-grady1);
subplot(2,2,3);
imshow(grady3,[]);
title('Substract Vertical');
 
gradxy3=abs(gradxy2-gradxy1);
subplot(2,2,4);
imshow(gradxy3,[]);
title('Substract Combination');
%=========================================================
% ---------------- question 3 ------------------------

% 1X2 Operator
No9 = figure;

%thresh_5 = compute_thresh_percent(gradxy1, 0.05);
%I3_thresh_5 = img_edge_thresh(gradxy1,thresh_5);

thresh_25 = compute_thresh_percent(gradxy1, 0.25);
I3_thresh_25 = img_edge_thresh(gradxy1,thresh_25);

subplot(1,2,1);
imshow(I3_thresh_25,[]);
title('1X2 Combination');

subplot(1,2,2);
histogram(I3_thresh_25);
title('1X2 Histogram');

% Sobel Operator
No10 = figure;

%thresh_5 = compute_thresh_percent(gradxy2, 0.05);
%I4_thresh_5 = img_edge_thresh(gradxy2,thresh_5);
thresh_25 = compute_thresh_percent(gradxy2, 0.25);
I4_thresh_25 = img_edge_thresh(gradxy2,thresh_25);

subplot(1,2,1);
imshow(I4_thresh_25,[]);
title('Sobel Combination');

subplot(1,2,2);
histogram(I4_thresh_25);
title('Sobel Histogram');


% Adaptive Sobel Threshold
T = adaptthresh(I2);
BW = imbinarize(I2,T);
No11 = figure();
imshow(BW);
title('Adaptive Threshold');

%=============================================================
% ---------------- qusestion 5------------------------

K1x2 = [-1 1];

K3x3 = [1 0 -1;2 0 -2;1 0 -1];

K5x5 = [2 1 0 -1 -2;3 2 0 -2 -3;4 3 0 -3 -4;3 2 0 -2 -3;2 1 0 -1 -2];

K7x7 = [3 2 1 0 -1 -2 -3;4 3 2 0 -2 -3 -4;5 4 3 0 -3 -4 -5;6 5 4 0 -4 -5 -6;5 4 3 0 -3 -4 -5;4 3 2 0 -2 -3 -4;3 2 1 0 -1 -2 -3];

No12 = figure;

G1x2 = edge_detection(I2, K1x2);
subplot(2,2,1);
imshow(G1x2,[]);
title("1x2")

G3x3 = edge_detection(I2, K3x3);
subplot(2,2,2);
imshow(G3x3,[]);
title("3x3")

G5x5 = edge_detection(I2, K5x5);
subplot(2,2,3);
imshow(G5x5,[]);
title("5x5")

G7x7 = edge_detection(I2, K7x7);
subplot(2,2,4);
imshow(G7x7,[]);
title("7x7")

% ---------------- Task 5 ------------------------

CR =uint8(I1(:,:,1));
CG =uint8(I1(:,:,2));
CB =uint8(I1(:,:,3));

Sobel = [-1 -2 -1; 0 0 0; 1 2 1];

No13 = figure();

subplot(2,2,1);
CRS = edge_detection(CR,Sobel);
imshow(CRS);
title('Red Sobel');

subplot(2,2,2);
CGS = edge_detection(CG,Sobel);
imshow(CGS);
title('Green Sobel');

subplot(2,2,3);
CBS = edge_detection(CB,Sobel);
imshow(CBS);
title('Blue Sobel');

subplot(2,2,4);
CS = 0.299 * CRS + 0.587 * CGS + 0.114 * CBS;
imshow(CS);
title('Combination Sobel');

No14 = figure;
subplot(2,2,1);
imshow(CS);
title('Edge map of combination');

subplot(2,2,2);
GCS = edge_detection(I2,Sobel);
imshow(GCS);
title('Edge map ofintensity image');

% RGB Combination subtract intensity image
subplot(2,2,3);
S = CS - GCS;
imshow(S);
title('RGB subtract Intensity')


No15 = figure;
subplot(2,2,1);
GCS = edge_detection(I2,Sobel);
imshow(GCS);
title('Edge map of intensity image');

subplot(2,2,2);
imshow(CS);
title('Edge map of combination');


% Color edge
subplot(2,2,3);
imshow(uint8(255 * cat(3, CRS, CGS, CBS)));
title('Color edge map')

No16 = figure;

thresholdR = imbinarize(CRS, adaptthresh(CRS));
thresholdG = imbinarize(CGS, adaptthresh(CGS));
thresholdB = imbinarize(CBS, adaptthresh(CBS));

CEM = uint8(255 * cat(3, thresholdR, thresholdG, thresholdB));

imshow(CEM);
title('Threshold color edge map')
%=========================================================
% ---------------- Functions ------------------------

function res = img_thresh(Img, thresh)
    res = Img;
    [ROWS, COLS, CHANNELS] = size(Img);
    for i=1 : ROWS
        for j=1 : COLS
            for m=1 : CHANNELS
              if Img(i,j,m) > thresh
                  res(i,j,m) = 255;
              else
                  res(i,j,m) = 0;
              end
            end
        end
    end
end

function thresh = compute_thresh_percent(Img, T)
    max_val = max(Img(:));
    thresh = max_val * T;
end

function res = img_edge_thresh(Img, thresh)
    res = Img;
    [ROWS, COLS, CHANNELS] = size(Img);
    for i=1 : ROWS
        for j=1 : COLS
            for m=1 : CHANNELS
              if Img(i,j,m) > thresh
                  res(i,j,m) = 1;
              else
                  res(i,j,m) = 0;
              end
            end
        end
    end
end

function res = edge_detection(img,kernel)
    Gx = kernel;
    Gy = Gx';
    tic
    gradx = norm(abs(conv2(img, Gx, 'same')));
    grady = norm(abs(conv2(img, Gy, 'same')));
    res = gradx+grady;
    toc
end

function normGrad = norm(grad)
    min_val = min(grad(:));
    max_val = max(grad(:));
    normGrad = (grad - min_val) ./ ( max_val - min_val);
end
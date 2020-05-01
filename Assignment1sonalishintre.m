%===================================================
% Computer Vision Programming Assignment 1
% @Zhigang Zhu, 
% City College of New Yor
%===================================================
% Sonali Shintre
% ID 7602

% ---------------- Step 1 ------------------------
% Read in an image, get information
% type help imread for more information
InputImage = 'http://ccvcl.org/wp-content/uploads/2019/09/IDPicture.png';
%InputImage = 'IDPicture.bmp'; 
%OutputImage1 = 'IDPicture_bw.bmp';

C1 = imread(InputImage);
[ROWS, COLS, CHANNELS] = size(C1);

% ---------------- Step 2 ------------------------
% If you want to display the three separate bands
% with the color image in one window, here is 
% what you need to do
% Basically you generate three "color" images
% using the three bands respectively
% and then use [] operator to concatenate the four images
% the orignal color, R band, G band and B band

% First, generate a blank image. Using "uinit8" will 
% give you an image of 8 bits for each pixel in each channel
% Since the Matlab will generate everything as double by default
CR1 =uint8(zeros(ROWS, COLS, CHANNELS));

% Note how to put the Red band of the color image C1 into 
% each band of the three-band grayscale image CR1
for band = 1 : CHANNELS
    CR1(:,:,band) = (C1(:,:,1));
end

% Do the same thing for G
CG1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS
    CG1(:,:,band) = (C1(:,:,2));
end

% and for B
CB1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS
    CB1(:,:,band) = (C1(:,:,3));
end

% Whenever you use figure, you generate a new figure window 
No1 = figure;  % Figure No. 1

%This is what I mean by concatenation
disimg = [C1, CR1;CG1, CB1]; 

% Then "image" will do the display for you!
image(disimg);
title('Orginal image')

% ---------------- Step 3 ------------------------
% Now we can calculate its intensity image from 
% the color image. Don't forget to use "uint8" to 
% covert the double results to unsigned 8-bit integers

I1    = uint8(round(sum(C1,3)/3));

% You can definitely display the black-white (grayscale)
% image directly without turn it into a three-band thing,
% which is a waste of memeory space

No2 = figure;  % Figure No. 2
image(I1);

% If you just stop your program here, you will see a 
% false color image since the system need a colormap to 
% display a 8-bit image  correctly. 
% The above display uses a default color map
% which is not correct. It is beautiful, though

% ---------------- Step 4 ------------------------
% So we need to generate a color map for the grayscale
% I think Matlab should have a function to do this,
% but I am going to do it myself anyway.

% Colormap is a 256 entry table, each index has three entries 
% indicating the three color components of the index

MAP =zeros(256, 3);

% For a gray scale C[i] = (i, i, i)
% But Matlab use color value from 0 to 1 
% so I scale 0-255 into 0-1 (and note 
% that I do not use "unit8" for MAP

for i = 1 : 256  % a comma means pause 
    for band = 1:CHANNELS
        MAP(i,band) = (i-1)/255;
    end 
end

%call colormap to enfore the MAP
colormap(MAP);

% I forgot to mention one thing: the index of Matlab starts from
% 1 instead 0.

% Is it correct this time? Remember the color table is 
% enforced for the current one, which is  the one we 
% just displayed.

% You can test if I am right by try to display the 
% intensity image again:

No3 = figure; % Figure No. 3
image(I1);


% See???
% You can actually check the color map using 
% the edit menu of each figure window

% ---------------- Step 5 ------------------------
% Use imwrite save any image
% check out image formats supported by Matlab
% by typing "help imwrite
% imwrite(I1, OutputImage1, 'BMP');


% ---------------- Step 6 and ... ------------------------
% Students need to do the rest of the jobs from c to g.
% Write code and comments - turn it in both in hard copies and 
% soft copies (electronically)

% ----------question 3-----------------------
%Intensity Image
%I = 0.299R + 0.587G + 0.114B
Intensity_image = 0.299*CR1 + 0.587*CG1 + 0.114*CB1;
No4 = figure;
image(Intensity_image);
title('Intensity Image');

% -----------question 4---------------------------
% Intensity image should have 256 gray levels.
% with K=4, 16, 32, 64    
% pixels whose values are below 128 are turned to 0,otherwise to 255. 
pixels= double(Intensity_image)/255;

%k=4
k4= uint8(pixels*4);
k4 = double(k4)/4;
No5 = figure();
image(k4);
title('k level = 4');

%k=16
k16= uint8(pixels*16);
k16 = double(k16)/16;
No6 = figure();
image(k16);
title('k level = 16');

%k=32
k32= uint8(pixels*32);
k32 = double(k32)/32;
No7 = figure();
image(k32);
title('k level = 32');

%k=64
k64= uint8(pixels*64);
k64 = double(k64)/64;
No8 = figure();
image(k64);
title('k level = 64');

% All images can repesented in an 2 dim array
No9 = figure;
subplot(2, 2, 1);
imshow(k4);
title('K = 4');

subplot(2, 2, 2);
imshow(k16);
title('K = 16');

subplot(2, 2, 3);
imshow(k32);
title('K = 32');

subplot(2, 2, 4);
imshow(k64);
title('K = 64');

%-------------question5---------------
% Now for the color image
% K level CK(x,y)= (R'(x,y), G'(x,y), B'(x,y)) 
%  K=2 and 4 each band 
colorpixels= double(C1)/255;

% k=2 color(kc)
kc2= uint8(colorpixels*2);
kc2 = double(kc2)/2;
No10 = figure();
image(kc2);
title('k level = 2');

%k=4
kc4= uint8(colorpixels*4);
kc4 = double (kc4)/4;
No11 = figure();
image(kc4);
title('k level = 4');

No12 = figure();
subplot(1, 2, 1);
imshow(kc2);
title('K = 2');
subplot(1, 2, 2);
imshow(kc4);
title('K = 4');
% we can conclude that as the k value increase image will be more clearly
% both in color(orginal image) and in the intensity image

%----------------question6---------------
%Quantize  the original three-band color image C1(x,y) 
% into a color image CL(x,y)= (R'(x,y), G'(x,y), B'(x,y)) (with a logarithmic function) , 
% the output range is still 0 – 255. Note that when I = 0, I' = 0 .
%I' =C ln (I+1) ( for each band), =C ln (I+1) ( for each band),


C= double (100)/log(100);
logimage= double(C1);
logimage= C * log(logimage +1);
logimage = uint8(logimage);

No14 = figure();
imshow(logimage);
title("log image");




 


%% ===========================================================
%% CSC I6716 Computer Vision
%% @ Zhigang Zhu, CCNY
%% Homework 4 - programming assignment:
%% Fundamental Matrix and Feature-Based Stereo Matching
%%
%% Name:Yaohua Chang
%% ID:23667288
%%
%% Note: Please do not delete the commented part of the code.
%% I am going to use it to test your program
%% =============================================================

%% Read in two images
imgl = imread('pic410.png');
imgr = imread('pic430.png');

%% display image pair side by side
[ROWS, COLS, CHANNELS] = size(imgl);
disimg = [imgl imgr];
image(disimg);

% You can change these two numbers,
% but use the variables; I will use them
% to load my data to test your algorithms

% Total Number of control points
Nc = 12;
% Total Number of test points
Nt = 4;

%% After several runs, you may want to save the point matches
%% in files (see below and then load them here, instead of
%% clicking many point matches every time you run the program

load pl.mat pl;
load pr.mat pr;

%% interface for picking up both the control points and
%% the test points

cnt = 1;
hold;

while(cnt <= Nc+Nt)

%% size of the rectangle to indicate point locations
dR = 50;
dC = 50;

%% pick up a point in the left image and display it with a rectangle....
%%% if you loaded the point matches, comment the point picking up (3 lines)%%%
% [X, Y] = ginput(1);
% Cl = X(1); Rl = Y(1);
% pl(cnt,:) = [Cl Rl 1];

%% and draw it
Cl= pl(cnt,1);  Rl=pl(cnt,2);
rectangle('Curvature', [0 0], 'Position', [Cl Rl dC dR]);

%% and then pick up the correspondence in the right image
%%% if you loaded the point matches, comment the point picking up (three lines)%%%

% [X, Y] = ginput(1);
% Cr = X(1); Rr = Y(1);
% pr(cnt,:) = [Cr-COLS Rr 1];

%% draw it
Cr=pr(cnt,1)+COLS; Rr=pr(cnt,2);
rectangle('Curvature', [0 0], 'Position', [Cr Rr dC dR]);
% plot(Cr+COLS,Rr,'r*');
drawnow;

cnt = cnt+1;
end

%% Student work (1a) NORMALIZATION: Page 156 of the textbook and Ex 7.6
%% --------------------------------------------------------------------
%% Normalize the coordinates of the corresponding points so that
%% the entries of A are of comparable size
%% You do not need to do this, but if you cannot get correct
%% result, you may want to use this
centroid = mean(pl,2);
dists = sqrt(sum((pl - repmat(centroid,1,size(pl,2))).^2,1));
scale = sqrt(2)/mean(dists);
norm_pl_T = scale * [1 0 -centroid(1);...
                  0 1 -centroid(2);...
                  0 0 1/scale];
norm_pl =  norm_pl_T * pl';
 
centroid = mean(pr,2);
dists = sqrt(sum((pr - repmat(centroid,1,size(pr,2))).^2,1));
scale = sqrt(2)/mean(dists);
norm_pr_T = scale * [1 0 -centroid(1);...
                  0 1 -centroid(2);...
                  0 0 1/scale];
norm_pr =  norm_pr_T * pr';

%% END NORMALIZATION %%

%% Student work: (1b) Implement EIGHT_POINT algorithm, page 156
%% --------------------------------------------------------------------
%% Generate the A matrix
A = [ repmat(norm_pr(1,:)',1,3) .* norm_pl', repmat(norm_pr(2,:)',1,3) .* norm_pl', norm_pl(1:3,:)'];


%% Singular value decomposition of A
[U,S,V] = svd(A);

%% the estimate of F
F_norm = reshape(V(:,end),3,3)';
[uf,sf,vf] = svd(F_norm);
F_norm_prime = uf*diag([sf(1) sf(5) 0])*(vf');

% Undo the coordinate normalization if you have done normalization
F= norm_pr_T' * F_norm_prime * norm_pl_T;


%% END of EIGHT_POINT

%% Draw the epipolar lines for both the controls points and the test
%% points, one by one; the current one (* in left and line in right) is in
%% red and the previous ones turn into blue

%% I suppose that your Fundamental matrix is F, a 3x3 matrix

%% Student work (1c): Check the accuray of the result by
%% measuring the distance between the estimated epipolar lines and
%% image points not used by the matrix estimation.
%% You can insert your code in the following for loop
error = zeros(Nc+Nt,1);
for cnt=1:1:Nc+Nt
  an = F*pl(cnt,:)';
  x = 0:COLS;
  y = -(an(1)*x+an(3))/an(2);

  hold on
  plot(pl(cnt,1),pl(cnt,2),'r*');
  plot(x,y,'Color', 'r');
  %[X, Y] = ginput(1); %% the location doesn't matter, press mouse to continue...
  
  an2 = F*pr(cnt,:)';
  y = -(an2(1)*x+an2(3))/an2(2);
  x = x + COLS;
  hold on
  plot(pr(cnt,1)+COLS,pr(cnt,2),'b*');
  line(x,y,'Color', 'b');

  error(cnt) = abs(an(2)*pr(cnt,2)+an(1)*pr(cnt,1)+an(3))/sqrt(an(1)^2+an(2)^2);
end

%% Save the corresponding points for later use... see discussions above
%save pr.mat pr;
%save pl.mat pl;

%% Save the F matrix in ascii
save F.txt F -ASCII

% Student work (1d): Find epipoles using the EPIPOLES_LOCATION algorithm page. 157
%% --------------------------------------------------------------------
[u,s,v] = svd(F);
el = v(:,end);
er = u(:,end);


%% save the eipoles

save eR.txt er -ASCII;
save eL.txt el -ASCII;

% Student work (2). Feature-based stereo matching
%% --------------------------------------------------------------------
%% Try to use the epipolar geometry derived from (1) in searching
%% correspondences along epipolar lines in Question (2). You may use
%% a similar interface  as I did for question (1). You may use the point
%% match searching algorithm in (1) (if you have done so), but this
%% time you need to constrain your search windows along the epipolar lines.

% Read the stereo images.
I1 = rgb2gray(imgl);
I2 = rgb2gray(imgr);

figure(2);
disimg = [imgl imgr];
imshow(disimg);
title('Select points on the first image');

% user interface that allows a user to select a point on the first image, say by a mouse click
cnt = 1;
hold;
while(cnt <= Nc+Nt)
    
% size of the rectangle to indicate point locations
dR = 50;
dC = 50;

% pick up a point in the left image and display it with a rectangle....
[X, Y] = ginput(1);
Cl = X(1);
Rl = Y(1);
pl(cnt,:) = [Cl Rl 1];

%and draw it
Cl = pl(cnt,1);
Rl = pl(cnt,2);
rectangle('Curvature', [0 0], 'Position', [Cl Rl dC dR]);
plot(Cl,Rl,'b*');
drawnow;

% crosshair points on right image
[X,Y] = crosshair_fun(Cl);
Cr = X(1);
Rr = Y(1);
pr(cnt,:) = [Cr Rr 1];

% [X, Y] = ginput(1);
% Cr = X(1); Rr = Y(1);
% pr(cnt,:) = [Cr-COLS Rr 1];

%% draw it
Cr = pr(cnt,1)+COLS;
Rr = pr(cnt,2);
rectangle('Curvature', [0 0], 'Position', [Cr Rr dC dR]);
plot(Cr,Rr,'r*');
drawnow;

cnt = cnt+1;
end

matchedPoints1 = pl(:, 1:2);
matchedPoints2 = pr(:, 1:2);

% Places the two images side by side
stackedImage = cat(2, imgl, imgr);
figure(3);
imshow(stackedImage);
width = size(imgl, 2);
hold on;
numPoints = size(matchedPoints1, 1);
% Note, we must offset by the width of the image
for i = 1 : numPoints
    plot(matchedPoints1(i, 1), matchedPoints1(i, 2), 'b*', matchedPoints2(i, 1) + width, matchedPoints2(i, 2), 'r*');
    line([matchedPoints1(i, 1) matchedPoints2(i, 1) + width], [matchedPoints1(i, 2) matchedPoints2(i, 2)], 'Color', 'yellow');
end
title('Automatically search point matches');

% crosshair points
function [X,Y] = crosshair_fun(x)
    X = abs(sin(x)) * x;
    Y = abs(cos(x)) * x;
end

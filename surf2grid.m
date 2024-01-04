function [Xq,Yq,Vq] = surf2grid(X,Y,V,stepSize,shpAlpha)
%Interpolates variably spaced surface mesh onto evenly space grid
%   X,Y = spacial locations of each variably mapped point
%   V = vector of +/- distance from 0 plane for each X,Y point
%   stepSize = distance spacing of new grid

% get upper and lower bounds of X,Y
min_X=floor(min(X(:))); min_Y=floor(min(Y(:)));
max_X=ceil(max(X(:))); max_Y=ceil(max(Y(:)));

[Xq,Yq] = meshgrid(min_X:stepSize:max_X, min_Y:stepSize:max_Y);

Vq = griddata(X,Y,V,Xq,Yq);

shp = alphaShape(X(:),Y(:),shpAlpha);
tf = inShape(shp,Xq,Yq);

Vq(~tf) = nan;

Xq = Xq(~all(isnan(Vq),2),:);
Yq = Yq(~all(isnan(Vq),2),:);
Vq = Vq(~all(isnan(Vq),2),:);

Xq = Xq(:,~all(isnan(Vq)));
Yq = Yq(:,~all(isnan(Vq)));
Vq = Vq(:,~all(isnan(Vq)));

end


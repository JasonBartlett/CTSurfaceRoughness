function [boxCounts,boxRatios] = irdbc(imageIn,boxRatios,graylevels,zerothLevel)
%% Calculates fractal dimension of grayscale images using
% Integer Ratio Differential Box Counting Method (IRDBC)
% Long, Min, and Fei Peng. A Box-Counting Method with Adaptable Box Height
% for Measuring the Fractal Feature of Images.
% Radioengineering, vol. 22, no. 1, 2013, pp. 208.
% imageIn: A MxN grayscale image matrix
% boxRatios: A vector of integers values whose range is from 2,...,n
% graylevels: Maximum number of possible graylevels eg 256, 128, etc
% zerothLevel: gray level equivelent to the zero surface deviation

% checks and makes sure imageIn is a matrix
if ~ismatrix(imageIn)
    error('Error. \nInput must be a MxN matrix')
end

% makes sure matrix is in double format for calculations
if ~isa(imageIn,'double')
    imageIn=double(imageIn);
end

% M=row count (Y dimension); N=column count (X dimension)
[M,N]=size(imageIn);

% identify total number of null cells on surface map image
imageNaNCount=sum(isnan(imageIn(:)));

% 2 <= r <= Q; Q=min(M,N,(MxN-NaN)^1/3);
% Q = min([floor((M*N)^(1/3)),M,N,floor(((M*N)-imageNaNCount)^(1/3))]);
Q = min([floor((M*N)^(1/3)),M,N]);

% trim boxRatios to only use values between 2 and Q inclusive
r=[boxRatios((boxRatios<=Q))];%,Q];

% create boxCounts vector (should be same size as r)
boxCounts=r.*0;

% determine variable boxsizes vectors (m,n) for used values of r
m=fix(M./r);
n=fix(N./r);

% calculate box area m and n for all values of r
mn=m.*n;

% box height p for for each r
p=graylevels./r;

% loop through system and calculate Nr(i,j)=sum(n(i,j))

for k=1:length(r)
    % set boxcount to 0
    Bc=0;

    %loop through image for each ratio size m,n
    for i=1:n(k):N
        for j=1:m(k):M

            % box selection
            if M>=j+m(k)-1 && N>=i+n(k)-1 % M=mr(k) & N=nr(k)
                currentBox=imageIn(j:j+m(k)-1,i:i+n(k)-1);
            elseif M<j+m(k)-1 && N>=i+n(k)-1 % M=mr(k) & N>nr(k)
                currentBox=imageIn(j:end,i:i+n(k)-1);
            elseif M>=j+m(k)-1 && N<i+n(k)-1 % M>mr(k) & N=nr(k)
                currentBox=imageIn(j:j+m(k)-1,i:end);
            elseif M<j+m(k)-1 && N<i+n(k)-1 % M>mr(k) & N>nr(k)
                currentBox=imageIn(j:end,i:end);
            end

            % area ratio n(i,j)/mn(k)
            nanCount=sum(isnan(currentBox(:)));

            % calculate actual area of coverage (boxSize-nan)
            Ar=((size(currentBox,1)*size(currentBox,2))-nanCount)/mn(k);

            % find minimum and maximum gray values (I) within selection
            maxI=max(currentBox(:),[],'omitnan');
            minI=min(currentBox(:),[],'omitnan');

            % see if range of minI:maxI cross the zeroth line
            % if maxI is below zeroth level set maxI to zeroth level value
            if maxI<zerothLevel || isnan(maxI)
                maxI=zerothLevel;
            end

            % if minI is above zeroth level set minI to zeroth level value
            if minI>zerothLevel || isnan(minI)
                minI=zerothLevel;
            end


            % calculate box grayscale intensity delta
            dI=abs(maxI-minI);

            % get fractional boxCount
            % sum Nr(i,j)
            Bc=Bc+(dI/p(k)+1)*Ar;

        end
    end

    boxCounts(k)= Bc;

end

boxRatios=r;

end
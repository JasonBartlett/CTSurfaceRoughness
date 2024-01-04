function [surfaceMeshData,subjectDataTable] = extractSurfaces(segmentData,subjectDataTable,voxelDimensions,varargin)
% converts original xyz segment slice coordinates into mesh surfaces
% surface meshes exported shape, curvature and total surface plots
% segmentData = struct containing airway tree segment surface data from XML extraction
% subjectDataTable = raw data table for surface stats
% voxelDimensions = CT voxel scales
% decimals = how many decimal points to use for rounding purposes
% alignment = 'fixed': used fixed points instead of top segment and bottom of segment

surfaceMeshData=struct;
decimalPlaces=4; % default decimal place
alignment='outerFirst'; % default alignment
surfaceRange=1:2; % process inner and outer surfaces

if ~isempty(varargin)
    if any(strcmp(varargin,'decimals'))
        decimalPlaces=varargin{find(strcmp(varargin,'decimals'))+1};
    end

    if any(strcmp(varargin,'alignment'))
        alignment=varargin{find(strcmp(varargin,'alignment'))+1}{:};
        if strcmp(alignment,'fixed')
            tracheaFixedRotationPoints=varargin{find(strcmp(varargin,'alignment'))+2};
        end
    end

end

% Segment|Link IDs
segIds=subjectDataTable.segId;
nPoints=73;

surfaceMeshData.subjectDataTable=subjectDataTable;

for kSeg=1:length(segIds)

    % check to see if segement data exisits
    if isfield(segmentData.(segIds{kSeg}),{'ConnectedPointXYZ'}) && ...
            isfield(segmentData.(segIds{kSeg}),{'ConnectedPointData'}) && ...
            sum(ismember(segmentData.(segIds{kSeg}...
            ).ConnectedPointData.Properties.VariableNames,'xC'))>0

        ConnectedPointXYZ=segmentData.(segIds{kSeg}).ConnectedPointXYZ;
        ConnectedPointData=segmentData.(segIds{kSeg}).ConnectedPointData;
        ConnectedPointData=ConnectedPointData(~isnan(ConnectedPointData.xC),:);

        % identify all connectedPoint Ids (aka Slices) {CPId_#)
        sliceIds=fieldnames(ConnectedPointXYZ);
        sliceCount=length(sliceIds);

        % determine which middle slices to use for Equivelent Circle
        numberOfAvgPoints=ceil(sliceCount/3);

        % get perimeter data for middle third to determine average diameter
        % (based on VIDAs Method)

        if mod((sliceCount-numberOfAvgPoints),2)==0

            middleThirdStart=(sliceCount-numberOfAvgPoints)/2+1;
            middleThirdEnd=sliceCount-(sliceCount-numberOfAvgPoints)/2;

        else

            middleThirdStart=floor((sliceCount-numberOfAvgPoints)/2)+1;
            middleThirdEnd=sliceCount-ceil((sliceCount-numberOfAvgPoints)/2);

        end

        % get perimeter data for each slice
        perimeterInner=ConnectedPointData.innerPerimeter(middleThirdStart:middleThirdEnd);
        perimeterOuter=ConnectedPointData.outerPerimeter(middleThirdStart:middleThirdEnd);

        % convert average of slices perimeters to equivelent circle diamter and radius
        eqCircleDiamInner=round(mean(perimeterInner)/(pi),decimalPlaces)';
        eqCircleDiamOuter=round(mean(perimeterOuter)/(pi),decimalPlaces)';

        eqCircleRadiusInner=round(mean(perimeterInner)/(2*pi),decimalPlaces)';
        eqCircleRadiusOuter=round(mean(perimeterOuter)/(2*pi),decimalPlaces)';

        % create coordinate arrays to be used for Point Clouds
        [innerCoordinates,outerCoordinates]=deal(ones(nPoints,3,length(sliceIds)));

        % identify mapped slices
        xyzCentresLUT=~isnan(ConnectedPointData.xC);

        % check that identified slice count matches actual slice count and
        % has at least 3 segmental slices
        if sum(xyzCentresLUT,1)==length(sliceIds) && length(sliceIds)>=3

            % populate inner and outer coordinate arrays
            for k=1:sliceCount

                % inner points (xyz)
                innerCoordinates(:,:,k)=[table2array(...
                    ConnectedPointXYZ.(sliceIds{k}).inner);...
                    table2array(... % wrap back around to first point closing the circle
                    ConnectedPointXYZ.(sliceIds{k}).inner(1,:))];

                % outer points (xyz)
                outerCoordinates(:,:,k)=[table2array(...
                    ConnectedPointXYZ.(sliceIds{k}).outer);...
                    table2array(... % wrap back around to first point closing the circle
                    ConnectedPointXYZ.(sliceIds{k}).outer(1,:))];

            end

            % entry point correction for trachea centerline starting point
            % to account for the slices which are cut off due to the
            % angular plane of the entry slices not being perfectly aligned
            % to the xy-planes. Ignores any slices not fully contained within
            % the positive z-plane when identifying the first center point
            % based on the selected surface.

            firstNonZeroInnerSlice=find(any(innerCoordinates(:,3,:)<0),1,'last')+1;
            firstNonZeroOuterSlice=find(any(outerCoordinates(:,3,:)<0),1,'last')+1;

            if isempty(firstNonZeroInnerSlice)
                firstNonZeroInnerSlice=1;
            end

            if isempty(firstNonZeroOuterSlice)
                firstNonZeroOuterSlice=1;
            end

            % reshape coordinates arrays into xyz planes whose size is 73 x sliceCount
            inner_x=reshape(innerCoordinates(:,1,:),[nPoints,sliceCount]);
            inner_y=reshape(innerCoordinates(:,2,:),[nPoints,sliceCount]);
            inner_z=reshape(innerCoordinates(:,3,:),[nPoints,sliceCount]);

            outer_x=reshape(outerCoordinates(:,1,:),[nPoints,sliceCount]);
            outer_y=reshape(outerCoordinates(:,2,:),[nPoints,sliceCount]);
            outer_z=reshape(outerCoordinates(:,3,:),[nPoints,sliceCount]);

            % extract segments slice centers
            segmentCentresRaw=[ConnectedPointData.xC,ConnectedPointData.yC,ConnectedPointData.zC];

            % scale segment coordinate arrays to mm dimensions and locations
            inner_x_scale=inner_x.*voxelDimensions.x;
            inner_y_scale=inner_y.*voxelDimensions.y;
            inner_z_scale=inner_z.*voxelDimensions.z;

            % scale segment coordinate arrays to mm dimensions and locations
            outer_x_scale=outer_x.*voxelDimensions.x;
            outer_y_scale=outer_y.*voxelDimensions.y;
            outer_z_scale=outer_z.*voxelDimensions.z;

            % recombine linearized scaled coordinate arrays
            innerSegmentScaled=[inner_x_scale(:),inner_y_scale(:),inner_z_scale(:)];
            outerSegmentScaled=[outer_x_scale(:),outer_y_scale(:),outer_z_scale(:)];

            % create unrotated & scaled segment PointCloud
            cmatrix = ones(size(outerSegmentScaled)).*[0 0 1];
            outerPCScaled=pointCloud(outerSegmentScaled,'color',cmatrix);
            innerPCScaled=pointCloud(innerSegmentScaled,'color',cmatrix);

            % scale segments slice centers
            segmentCentresScaled=segmentCentresRaw.*([voxelDimensions.x,voxelDimensions.y,voxelDimensions.z]);

            % Set slices to use for segment alignment

            if strcmp(subjectDataTable.anatomicalName(strcmp(subjectDataTable.segId,segIds{kSeg})),'Trachea')

                if strcmp(alignment,'outerFirst')

                    segAlignStart=firstNonZeroOuterSlice;
                    segAlignEnd=sliceCount;

                elseif strcmp(alignment,'innerFirst')

                    segAlignStart=firstNonZeroInnerSlice;
                    segAlignEnd=sliceCount;

                elseif strcmp(alignment,'middleThird')

                    segAlignStart=middleThirdStart;
                    segAlignEnd=middleThirdEnd;

                elseif strcmp(alignment,'fixed')

                    sliceIndex=floor(segmentCentresRaw(:,3))>=tracheaFixedRotationPoints(1)&floor(segmentCentresRaw(:,3))<=tracheaFixedRotationPoints(2);

                    segAlignStart=find(sliceIndex,1,'first') ;
                    segAlignEnd=find(sliceIndex,1,'last');

                    if segAlignStart<firstNonZeroOuterSlice
                        segAlignStart=firstNonZeroOuterSlice;
                    end

                elseif strcmp(alignment,'topbot')

                    segAlignStart=1;
                    segAlignEnd=sliceCount;

                end

            else

                if strcmp(alignment,'middleThird')

                    segAlignStart=middleThirdStart;
                    segAlignEnd=middleThirdEnd;

                else

                    % unless otherwise defined for any other segment
                    % align the segment to the the top and bottom
                    % center points of the segment

                    segAlignStart=1;
                    segAlignEnd=sliceCount;

                end
            end


            % extract center line for rigid transformation alignment
            centerLine2Align=segmentCentresScaled([segAlignStart,segAlignEnd],:)';
            lineCenters=sum(centerLine2Align,2)/2;

            % create vertical center line to align center line to
            centerLineFixed=centerLine2Align;
            centerLineFixed(1:2,:)=ones(2,2).*lineCenters(1:2);

            % create rigid transformation matrix for segment rotation vertical alignment
            centerLine2AlignPoints=[centerLine2Align(:,1),lineCenters,centerLine2Align(:,2)]';
            centerLineFixedPoints=[centerLineFixed(:,1),lineCenters,centerLineFixed(:,2)]';

            centerAlignmentCloudMoving=pointCloud(centerLine2AlignPoints);
            centerAlignmentCloudFixed=pointCloud(centerLineFixedPoints);

            tform = pcregistericp(centerAlignmentCloudMoving,centerAlignmentCloudFixed);

            % create centerline pointcloud and apply rotation matix
            centerAlignmentCloudMoving=pointCloud(centerLine2Align');
            lineRotated=pctransform(centerAlignmentCloudMoving,tform);
            centerLineCloud=pointCloud(segmentCentresScaled);
            centerLineRotated=pctransform(centerLineCloud,tform);

            % apply rotation matix to scaled segment pointcloud
            outerPCRotated=pctransform(outerPCScaled,tform);
            innerPCRotated=pctransform(innerPCScaled,tform);

            % extract segment and centerline pointcloud coordinates
            segmentFixedCenter=lineRotated.Location;
            centerLinesRotated=centerLineRotated.Location;

            innerRotatedPC=innerPCRotated.Location;
            outerRotatedPC=outerPCRotated.Location;

            % align center line coordinates to x=0 and y=0 and set very
            % small <10^-6 differences to 0
            centerLinesRotated(:,1)=centerLinesRotated(:,1)-segmentFixedCenter(1,1);
            centerLinesRotated(abs(centerLinesRotated(:,1))<1E-6,1)=0;

            centerLinesRotated(:,2)=centerLinesRotated(:,2)-segmentFixedCenter(1,2);
            centerLinesRotated(abs(centerLinesRotated(:,2))<1E-6,2)=0;

            zero_z=min(min(innerRotatedPC(:,3),[],'all'),min(outerRotatedPC(:,3),[],'all'));
            centerLinesRotated(:,3)=centerLinesRotated(:,3)-zero_z;
            centerLinesRotated(abs(centerLinesRotated(:,3))<1E-6,3)=0;

            % align center of segment to x=0
            innerRotatedPC(:,1)=innerRotatedPC(:,1)-segmentFixedCenter(1,1);
            % align center of segment to y=0
            innerRotatedPC(:,2)=innerRotatedPC(:,2)-segmentFixedCenter(1,2);
            % align bottom of segment to z=0
            innerRotatedPC(:,3)=innerRotatedPC(:,3)-zero_z;

            % align center of segment to x=0
            outerRotatedPC(:,1)=outerRotatedPC(:,1)-segmentFixedCenter(1,1);
            % align center of segment to y=0
            outerRotatedPC(:,2)=outerRotatedPC(:,2)-segmentFixedCenter(1,2);
            % align bottom of segment to z=0
            outerRotatedPC(:,3)=outerRotatedPC(:,3)-zero_z;

            % calculate four quadrant inverse angle in radians (atan(Y/X))
            innerSlicedRotatedSegments_x=reshape(innerRotatedPC(:,1),[73,sliceCount]);
            innerSlicedRotatedSegments_y=reshape(innerRotatedPC(:,2),[73,sliceCount]);
            innerSlicedRotatedSegments_z=reshape(innerRotatedPC(:,3),[73,sliceCount]);

            outerSlicedRotatedSegments_x=reshape(outerRotatedPC(:,1),[73,sliceCount]);
            outerSlicedRotatedSegments_y=reshape(outerRotatedPC(:,2),[73,sliceCount]);
            outerSlicedRotatedSegments_z=reshape(outerRotatedPC(:,3),[73,sliceCount]);

            innerSlicedRotatedSegments_x_zeroed=innerSlicedRotatedSegments_x-centerLinesRotated(:,1)';
            innerSlicedRotatedSegments_y_zeroed=innerSlicedRotatedSegments_y-centerLinesRotated(:,2)';

            outerSlicedRotatedSegments_x_zeroed=outerSlicedRotatedSegments_x-centerLinesRotated(:,1)';
            outerSlicedRotatedSegments_y_zeroed=outerSlicedRotatedSegments_y-centerLinesRotated(:,2)';

            innerThetaRad=atan2(innerSlicedRotatedSegments_x_zeroed,innerSlicedRotatedSegments_y_zeroed);
            outerThetaRad=atan2(outerSlicedRotatedSegments_x_zeroed,outerSlicedRotatedSegments_y_zeroed);

            % extract surface deviation topological values for
            % Shape, Curvature, Shape+Curvature
            % based on Rotated and scaled Segments

            % generate surface shape distortion by comparing each slices
            % shape to that of a perfectly round circle
            % by transforming each slice point such that the center of the
            % slice becomes x,y=0,0

            innerShapeDeviationsSlices=sqrt(innerSlicedRotatedSegments_x_zeroed(:,:,1).^2+...
                innerSlicedRotatedSegments_y_zeroed.^2)-...
                eqCircleRadiusInner;

            outerShapeDeviationsSlices=sqrt(outerSlicedRotatedSegments_x_zeroed.^2+...
                outerSlicedRotatedSegments_y_zeroed.^2)-...
                eqCircleRadiusOuter;

            % generate perfect cylindar point cloud for export and surface
            % curvature distortion mapping

            innerSurfacePerfect_x=sin(innerThetaRad).*eqCircleRadiusInner;
            innerSurfacePerfect_y=cos(innerThetaRad).*eqCircleRadiusInner;

            outerSurfacePerfect_x=sin(innerThetaRad).*eqCircleRadiusOuter;
            outerSurfacePerfect_y=cos(innerThetaRad).*eqCircleRadiusOuter;

            % generate surface curvature devations using a perfect cylinder
            % by applying the centerline distortion to the perfect cylinder
            % segment slices centerline points and recording the relative
            % distortion of each surface point

            innerSurfacePerfectCurved_x=innerSurfacePerfect_x+centerLinesRotated(:,1)';
            innerSurfacePerfectCurved_y=innerSurfacePerfect_y+centerLinesRotated(:,2)';

            outerSurfacePerfectCurved_x=outerSurfacePerfect_x+centerLinesRotated(:,1)';
            outerSurfacePerfectCurved_y=outerSurfacePerfect_y+centerLinesRotated(:,2)';

            % determine new vectors by moving perfect surface point to
            % x1,y1=(0,0) and ajusting x0,y0=(0-x1,0-y1) and matching
            % surface points of perfect curved surface x2,y2=(x2-x1,y2-y1)

            inner_zeroedVec_x_perf=0-innerSurfacePerfect_x;
            inner_zeroedVec_y_perf=0-innerSurfacePerfect_y;

            inner_zeroedVec_x_perfcurv=innerSurfacePerfectCurved_x-innerSurfacePerfect_x;
            inner_zeroedVec_y_perfcurv=innerSurfacePerfectCurved_y-innerSurfacePerfect_y;

            outer_zeroedVec_x_perf=0-outerSurfacePerfect_x;
            outer_zeroedVec_y_perf=0-outerSurfacePerfect_y;

            outer_zeroedVec_x_perfcurv=outerSurfacePerfectCurved_x-outerSurfacePerfect_x;
            outer_zeroedVec_y_perfcurv=outerSurfacePerfectCurved_y-outerSurfacePerfect_y;

            % calculate distance between perfect and perfect curved surface
            % points

            inner_hypotenuse_perfcurv=sqrt(inner_zeroedVec_x_perfcurv.^2 + inner_zeroedVec_y_perfcurv.^2);
            outer_hypotenuse_perfcurv=sqrt(outer_zeroedVec_x_perfcurv.^2 + outer_zeroedVec_y_perfcurv.^2);

            % determine angle between x0,y0 and x2,y2

            innerCurvatureAnglesRad = atan2(...
                inner_zeroedVec_x_perf.*inner_zeroedVec_y_perfcurv ...
                - inner_zeroedVec_y_perf.*inner_zeroedVec_x_perfcurv, ...
                inner_zeroedVec_x_perf.*inner_zeroedVec_x_perfcurv ...
                + inner_zeroedVec_y_perf.*inner_zeroedVec_y_perfcurv);

            outerCurvatureAnglesRad = atan2(...
                outer_zeroedVec_x_perf.*outer_zeroedVec_y_perfcurv ...
                - outer_zeroedVec_y_perf.*outer_zeroedVec_x_perfcurv, ...
                outer_zeroedVec_x_perf.*outer_zeroedVec_x_perfcurv ...
                + outer_zeroedVec_y_perf.*outer_zeroedVec_y_perfcurv);

            % determine amount that perfect curved surface devaiates
            % in the delta direction
            % (eg vertically from unrolled perfectly smooth surface)

            innerCurveDeviationsSlices=(cos(innerCurvatureAnglesRad(:,:)).*inner_hypotenuse_perfcurv(:,:)).*-1;
            outerCurveDeviationsSlices=(cos(outerCurvatureAnglesRad(:,:)).*outer_hypotenuse_perfcurv(:,:)).*-1;

            % generate actual surface curvature + distortions devations
            % using the perfect cylinder recording the relative
            % distortion of each surface point to the perfect cylinder

            % determine new vectors by moving perfect surface point to
            % x1,y1=(0,0) and ajusting x0,y0=(0-x1,0-y1) and matching
            % surface points of actual curved surface x2,y2=(x2-x1,y2-y1)

            inner_zeroedVec_x_surfcurv=innerSlicedRotatedSegments_x-innerSurfacePerfect_x;
            inner_zeroedVec_y_surfcurv=innerSlicedRotatedSegments_y-innerSurfacePerfect_y;

            outer_zeroedVec_x_surfcurv=outerSlicedRotatedSegments_x-outerSurfacePerfect_x;
            outer_zeroedVec_y_surfcurv=outerSlicedRotatedSegments_y-outerSurfacePerfect_y;

            % calculate distance between perfect and surface with curvature
            % equilvent points

            inner_hypotenuse_surfcurv=sqrt(inner_zeroedVec_x_surfcurv.^2 + inner_zeroedVec_y_surfcurv.^2);
            outer_hypotenuse_surfcurv=sqrt(outer_zeroedVec_x_surfcurv.^2 + outer_zeroedVec_y_surfcurv.^2);

            % determine angle between x0,y0 and x2,y2

            innerSurfCurvAnglesRad = atan2(...
                inner_zeroedVec_x_perf.*inner_zeroedVec_y_surfcurv ...
                - inner_zeroedVec_y_perf.*inner_zeroedVec_x_surfcurv, ...
                inner_zeroedVec_x_perf.*inner_zeroedVec_x_surfcurv ...
                + inner_zeroedVec_y_perf.*inner_zeroedVec_y_surfcurv);

            outerSurfCurvAnglesRad = atan2(...
                outer_zeroedVec_x_perf.*outer_zeroedVec_y_surfcurv ...
                - outer_zeroedVec_y_perf.*outer_zeroedVec_x_surfcurv, ...
                outer_zeroedVec_x_perf.*outer_zeroedVec_x_surfcurv ...
                + outer_zeroedVec_y_perf.*outer_zeroedVec_y_surfcurv);

            % determine amount that actual surface with curvature devaiates
            % in the delta direction
            % (eg vertically from unrolled perfectly smooth surface)

            innerTotalDeviationsSlices=(cos(innerSurfCurvAnglesRad(:,:)).*inner_hypotenuse_surfcurv(:,:)).*-1;
            outerTotalDeviationsSlices=(cos(outerSurfCurvAnglesRad(:,:)).*outer_hypotenuse_surfcurv(:,:)).*-1;


            % determine variable arc spacing for surface points
            % using chord length to arc length conversion

            innerArcLengthsSlices=innerThetaRad.*eqCircleRadiusInner;
            outerArcLengthsSlices=outerThetaRad.*eqCircleRadiusOuter;

            % square differences in distances between unrotated segment slices x,y,z points
            inner_x_diff_scale2=diff(inner_x_scale,1).^2;
            inner_y_diff_scale2=diff(inner_y_scale,1).^2;
            inner_z_diff_scale2=diff(inner_z_scale,1).^2;

            outer_x_diff_scale2=diff(outer_x_scale,1).^2;
            outer_y_diff_scale2=diff(outer_y_scale,1).^2;
            outer_z_diff_scale2=diff(outer_z_scale,1).^2;

            % determine euclidian distance between each of the slices 72 points
            innerDistBetweenPoints=sqrt(inner_x_diff_scale2+inner_y_diff_scale2+inner_z_diff_scale2);
            outerDistBetweenPoints=sqrt(outer_x_diff_scale2+outer_y_diff_scale2+outer_z_diff_scale2);

            % determine z distance between each of the slices 72 points
            innerDistBetweenZSlice=diff(innerSlicedRotatedSegments_z,1);
            outerDistBetweenZSlice=diff(outerSlicedRotatedSegments_z,1);

            % determine chordLenght using c^2-a^2=b^2 (hyp^2-opp^2=adj^2)
            innerChordLength=sqrt(innerDistBetweenPoints.^2-innerDistBetweenZSlice.^2);
            outerChordLength=sqrt(outerDistBetweenPoints.^2-outerDistBetweenZSlice.^2);

            % determine arc length distance between the each of the 72 points
            innerArcLengthVariable=asin(innerChordLength./(2*eqCircleRadiusInner)).*(2*eqCircleRadiusInner);
            outerArcLengthVariable=asin(outerChordLength./(2*eqCircleRadiusOuter)).*(2*eqCircleRadiusOuter);

            % sort baseline arc lengths and remove wrap around point
            [innerArcLengthSorted,innerSortLut]=sort(innerArcLengthsSlices(1:end-1,:),1,'descend');
            [outerArcLengthSorted,outerSortLut]=sort(outerArcLengthsSlices(1:end-1,:),1,'descend');

            % preallocate thetaRad, arcLengthVariable, radiiDeviationsSlices, zSlices for sorting

            innerThetaSorted=innerThetaRad(1:end-1,:);
            outerThetaSorted=outerThetaRad(1:end-1,:);

            innerSlicedRotatedSegments_x_Sorted=innerSlicedRotatedSegments_x(1:end-1,:);
            innerSlicedRotatedSegments_y_Sorted=innerSlicedRotatedSegments_y(1:end-1,:);

            innerSlicedPerfectSegments_x_Sorted=innerSurfacePerfect_x(1:end-1,:);
            innerSlicedPerfectSegments_y_Sorted=innerSurfacePerfect_y(1:end-1,:);

            innerArcLengthVariableSorted=innerArcLengthVariable;
            innerDeviationsSurfOnlySorted=innerShapeDeviationsSlices(1:end-1,:);
            innerDeviationsCurvOnlySorted=innerCurveDeviationsSlices(1:end-1,:);
            innerDeviationsSurfCurvSorted=innerTotalDeviationsSlices(1:end-1,:);
            inner_zSlicesSorted=innerSlicedRotatedSegments_z(1:end-1,:);

            outerSlicedRotatedSegments_x_Sorted=outerSlicedRotatedSegments_x(1:end-1,:);
            outerSlicedRotatedSegments_y_Sorted=outerSlicedRotatedSegments_y(1:end-1,:);

            outerSlicedPerfectSegments_x_Sorted=outerSurfacePerfect_x(1:end-1,:);
            outerSlicedPerfectSegments_y_Sorted=outerSurfacePerfect_y(1:end-1,:);

            outerArcLengthVariableSorted=outerArcLengthVariable;
            outerDeviationsSurfOnlySorted=outerShapeDeviationsSlices(1:end-1,:);
            outerDeviationsCurvOnlySorted=outerCurveDeviationsSlices(1:end-1,:);
            outerDeviationsSurfCurvSorted=outerTotalDeviationsSlices(1:end-1,:);
            outer_zSlicesSorted=outerSlicedRotatedSegments_z(1:end-1,:);

            % sort matrices according to sorted baseline arc lengths
            for kcol = 1:sliceCount

                % inner rotated segments
                innerSlicedRotatedSegments_x_Sorted(:,kcol)=...
                    innerSlicedRotatedSegments_x_Sorted(innerSortLut(:,kcol),kcol);
                innerSlicedRotatedSegments_y_Sorted(:,kcol)=...
                    innerSlicedRotatedSegments_y_Sorted(innerSortLut(:,kcol),kcol);

                % outer rotated segments
                outerSlicedRotatedSegments_x_Sorted(:,kcol)=...
                    outerSlicedRotatedSegments_x_Sorted(outerSortLut(:,kcol),kcol);
                outerSlicedRotatedSegments_y_Sorted(:,kcol)=...
                    outerSlicedRotatedSegments_y_Sorted(outerSortLut(:,kcol),kcol);


                % inner perfect segments
                innerSlicedPerfectSegments_x_Sorted(:,kcol)=...
                    innerSlicedPerfectSegments_x_Sorted(innerSortLut(:,kcol),kcol);
                innerSlicedPerfectSegments_y_Sorted(:,kcol)=...
                    innerSlicedPerfectSegments_y_Sorted(innerSortLut(:,kcol),kcol);

                % outer perfect segments
                outerSlicedPerfectSegments_x_Sorted(:,kcol)=...
                    outerSlicedPerfectSegments_x_Sorted(outerSortLut(:,kcol),kcol);
                outerSlicedPerfectSegments_y_Sorted(:,kcol)=...
                    outerSlicedPerfectSegments_y_Sorted(outerSortLut(:,kcol),kcol);                                

                % thetaRad

                innerThetaSorted(:,kcol)=innerThetaSorted(innerSortLut(:,kcol),kcol);
                outerThetaSorted(:,kcol)=outerThetaSorted(innerSortLut(:,kcol),kcol);

                % q1
                inner_zSlicesSorted(:,kcol) = inner_zSlicesSorted(innerSortLut(:,kcol),kcol);
                outer_zSlicesSorted(:,kcol) = outer_zSlicesSorted(outerSortLut(:,kcol),kcol);

                % q2
                innerArcLengthVariableSorted(:,kcol) = innerArcLengthVariable(innerSortLut(:,kcol),kcol);
                outerArcLengthVariableSorted(:,kcol) = outerArcLengthVariable(outerSortLut(:,kcol),kcol);

                % % deltas
                innerDeviationsSurfOnlySorted(:,kcol) = innerDeviationsSurfOnlySorted(innerSortLut(:,kcol),kcol);
                outerDeviationsSurfOnlySorted(:,kcol) = outerDeviationsSurfOnlySorted(outerSortLut(:,kcol),kcol);
                %
                innerDeviationsCurvOnlySorted(:,kcol) = innerDeviationsCurvOnlySorted(innerSortLut(:,kcol),kcol);
                outerDeviationsCurvOnlySorted(:,kcol) = outerDeviationsCurvOnlySorted(outerSortLut(:,kcol),kcol);
                %
                innerDeviationsSurfCurvSorted(:,kcol) = innerDeviationsSurfCurvSorted(innerSortLut(:,kcol),kcol);
                outerDeviationsSurfCurvSorted(:,kcol) = outerDeviationsSurfCurvSorted(outerSortLut(:,kcol),kcol);

            end

            % idenify smallest magnitude arc length location for all slices
            [innerArcMin,innerArcMinLoc]=min(abs(innerArcLengthSorted),[],1);
            [outerArcMin,outerArcMinLoc]=min(abs(outerArcLengthSorted),[],1);

            % Preallocate arcLengthVariableFinal memory space
            innerArcLengthVariableSortedFinal=innerArcLengthVariableSorted;
            outerArcLengthVariableSortedFinal=outerArcLengthVariableSorted;

            % recalculate radii arc spacings based on variable arc lengths distances
            for karc = 1:sliceCount

                innerSliceMinValue=innerArcMin(karc);
                innerSliceMinLoc=innerArcMinLoc(karc);
                % cumulatively subtract arc distances between points (0,-pi) from smallest magnitude value
                innerArcLengthVariableSortedFinal(1:innerSliceMinLoc-1,karc)...
                    =innerSliceMinValue-cumsum(innerArcLengthVariableSorted(1:innerSliceMinLoc-1,karc),'reverse');
                % replace middle point with smallest magnitude value
                innerArcLengthVariableSortedFinal(innerSliceMinLoc,karc)=innerSliceMinValue;
                % cumulatively add arc distances between points (0,pi) to smallest magnitude value
                innerArcLengthVariableSortedFinal(innerSliceMinLoc+1:end,karc)...
                    =innerSliceMinValue+cumsum(innerArcLengthVariableSorted(innerSliceMinLoc+1:end,karc));

                outerSliceMinValue=outerArcMin(karc);
                outerSliceMinLoc=outerArcMinLoc(karc);
                % cumulatively subtract arc distances between points (0,-pi) from smallest magnitude value
                outerArcLengthVariableSortedFinal(1:outerSliceMinLoc-1,karc)...
                    =outerSliceMinValue-cumsum(outerArcLengthVariableSorted(1:outerSliceMinLoc-1,karc),'reverse');
                % replace middle point with smallest magnitude value
                outerArcLengthVariableSortedFinal(outerSliceMinLoc,karc)=outerSliceMinValue;
                % cumulatively add arc distances between points (0,pi) to smallest magnitude value
                outerArcLengthVariableSortedFinal(outerSliceMinLoc+1:end,karc)...
                    =outerSliceMinValue+cumsum(outerArcLengthVariableSorted(outerSliceMinLoc+1:end,karc));

            end

            innerArcLengthVariableSortedFinal=flip(innerArcLengthVariableSortedFinal);
            outerArcLengthVariableSortedFinal=flip(outerArcLengthVariableSortedFinal);

            % compile segment point cloud data for export
            % Original = unrotated unscaled unaligned voxel scaled segment
            % Rotated = scaled rotated centered mm scaled segment

            % inner segments

            innerSegmentSurfaceOriginal=(zeros(72,sliceCount,3));
            innerSegmentSurfaceRotated=innerSegmentSurfaceOriginal;
            innerSegmentPerfect=innerSegmentSurfaceOriginal;

            innerSegmentSurfaceOriginal(:,:,1)=inner_x(1:72,:);
            innerSegmentSurfaceOriginal(:,:,2)=inner_y(1:72,:);
            innerSegmentSurfaceOriginal(:,:,3)=inner_z(1:72,:);

            innerSegmentSurfaceRotated(:,:,1)=innerSlicedRotatedSegments_x_Sorted;
            innerSegmentSurfaceRotated(:,:,2)=innerSlicedRotatedSegments_y_Sorted;
            innerSegmentSurfaceRotated(:,:,3)=inner_zSlicesSorted;

            innerSegmentPerfect(:,:,1)=innerSlicedPerfectSegments_x_Sorted;
            innerSegmentPerfect(:,:,2)=innerSlicedPerfectSegments_y_Sorted;
            innerSegmentPerfect(:,:,3)=inner_zSlicesSorted;

            % outer segments

            outerSegmentSurfaceOriginal=(zeros(72,sliceCount,3));
            outerSegmentSurfaceRotated=outerSegmentSurfaceOriginal;
            outerSegmentPerfect=outerSegmentSurfaceOriginal;

            outerSegmentSurfaceOriginal(:,:,1)=outer_x(1:72,:);
            outerSegmentSurfaceOriginal(:,:,2)=outer_y(1:72,:);
            outerSegmentSurfaceOriginal(:,:,3)=outer_z(1:72,:);

            outerSegmentSurfaceRotated(:,:,1)=outerSlicedRotatedSegments_x_Sorted;
            outerSegmentSurfaceRotated(:,:,2)=outerSlicedRotatedSegments_y_Sorted;
            outerSegmentSurfaceRotated(:,:,3)=outer_zSlicesSorted;

            outerSegmentPerfect(:,:,1)=outerSlicedPerfectSegments_x_Sorted;
            outerSegmentPerfect(:,:,2)=outerSlicedPerfectSegments_y_Sorted;
            outerSegmentPerfect(:,:,3)=inner_zSlicesSorted;
            
            % prepare surface plots for export
            innerShapeMesh=(zeros(72,sliceCount,3));
            outerShapeMesh=(zeros(72,sliceCount,3));

            % Basic mapping is ...

            % q1=z -> location of coordinate point along z-axis
            innerShapeMesh(:,:,1)=inner_zSlicesSorted;
            outerShapeMesh(:,:,1)=outer_zSlicesSorted;

            innerCurveMesh=innerShapeMesh;
            outerCurveMesh=outerShapeMesh;

            % q2=arcLengths -> arc dist location of coordinate point along perfect clyindar
            innerShapeMesh(:,:,2)=innerArcLengthVariableSortedFinal;
            outerShapeMesh(:,:,2)=outerArcLengthVariableSortedFinal;

            innerCurveMesh(:,:,2)=innerArcLengthSorted;
            outerCurveMesh(:,:,2)=outerArcLengthSorted;

            innerTotalMesh=innerShapeMesh;
            outerTotalMesh=outerShapeMesh;           

            % delta -> deviation of coordinate point from perfect clyindar
            innerShapeMesh(:,:,3)=innerDeviationsSurfOnlySorted;
            outerShapeMesh(:,:,3)=outerDeviationsSurfOnlySorted;

            innerCurveMesh(:,:,3)=innerDeviationsCurvOnlySorted;
            outerCurveMesh(:,:,3)=outerDeviationsCurvOnlySorted;

            innerTotalMesh(:,:,3)=innerDeviationsSurfCurvSorted;
            outerTotalMesh(:,:,3)=outerDeviationsSurfCurvSorted;

            % store exports into

            surfaceMeshData.subjectDataTable.sliceCount(kSeg)=sliceCount;
            surfaceMeshData.subjectDataTable.innerNonZeroSlice(kSeg)=firstNonZeroInnerSlice;
            surfaceMeshData.subjectDataTable.outerNonZeroSlice(kSeg)=firstNonZeroOuterSlice;
            surfaceMeshData.subjectDataTable.middleThirdStart(kSeg)=middleThirdStart;
            surfaceMeshData.subjectDataTable.middleThirdEnd(kSeg)=middleThirdEnd;
            surfaceMeshData.subjectDataTable.numberOfAvgPoints(kSeg)=numberOfAvgPoints;
            surfaceMeshData.subjectDataTable.innerEqCDiameter(kSeg)=round(eqCircleDiamInner,decimalPlaces);
            surfaceMeshData.subjectDataTable.innerEqCRadius(kSeg)=round(eqCircleRadiusInner,decimalPlaces);
            surfaceMeshData.subjectDataTable.outerEqCDiameter(kSeg)=round(eqCircleDiamOuter,decimalPlaces);
            surfaceMeshData.subjectDataTable.outerEqCRadius(kSeg)=round(eqCircleRadiusOuter,decimalPlaces);
            surfaceMeshData.subjectDataTable.alignment{kSeg}=alignment;
            surfaceMeshData.subjectDataTable.segAlignStart(kSeg)=segAlignStart;
            surfaceMeshData.subjectDataTable.segAlignEnd(kSeg)=segAlignEnd;
            surfaceMeshData.subjectDataTable.segTranslatex(kSeg)=round(segmentFixedCenter(1,1),decimalPlaces);
            surfaceMeshData.subjectDataTable.segTranslatey(kSeg)=round(segmentFixedCenter(1,2),decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).segmentRotationMatrix=tform;

            surfaceMeshData.(segIds{kSeg}).centerLinePointsOriginal=round(segmentCentresRaw,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).centerLinePointsScaled=round(centerLinesRotated,decimalPlaces);
            
            surfaceMeshData.(segIds{kSeg}).inner.CylThetaRad=round(innerThetaSorted,decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).inner.OriginalSegment=round(innerSegmentSurfaceOriginal,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).inner.ScaledSegment=round(innerSegmentSurfaceRotated,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).inner.PerfectSegment=round(innerSegmentPerfect,decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).inner.ShapeMesh=round(innerShapeMesh,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).inner.CurveMesh=round(innerCurveMesh,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).inner.TotalMesh=round(innerTotalMesh,decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).outer.CylThetaRad=round(outerThetaSorted,decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).outer.OriginalSegment=round(outerSegmentSurfaceOriginal,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).outer.ScaledSegment=round(outerSegmentSurfaceRotated,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).outer.PerfectSegment=round(outerSegmentPerfect,decimalPlaces);

            surfaceMeshData.(segIds{kSeg}).outer.ShapeMesh=round(outerShapeMesh,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).outer.CurveMesh=round(outerCurveMesh,decimalPlaces);
            surfaceMeshData.(segIds{kSeg}).outer.TotalMesh=round(outerTotalMesh,decimalPlaces);

        else
            if length(sliceIds)<3
                disp(strcat("  Error: less than 3 slices ",segIds{kSeg}))
            else
                disp(strcat("  Error: slice count mismatch ",segIds{kSeg}))
            end
        end

    else

        disp(strcat("  Error: segment data does not exist ",segIds{kSeg}))

    end

    disp(strcat(" - Finished processing segment: ",segIds{kSeg}))

end

disp(strcat("Finished extracting all surfaces"))


end

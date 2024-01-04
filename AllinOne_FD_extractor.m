clear; close; clc;

%% all in one pass for extracting trachea surfaces from raw VIDA output files
% requires VIDA XML and masks for aortic arch and lungs to use as anchor
% points.

% requires the following functions
% readairwayXML.m - converted xml to surface segment point clouds
% extractSurfaces.m - converts segment point clouds to surface mesh plots
% surf2grid.m - coverts surface mesh plots to 2d grayscale surface height maps
% irdbcRatioSelect.m - identifys ratio values to use based on maximizing box coverage
% irdbc.m - extracts Fractal dimensions from 2d grayscale surface height maps
% checkFolder.m - checks for and creates empty folders if missing for final output

%% Set Variables
importRootFolder='R:\kirby_group\CanCOLD\apks'; % master subject directory
exportXmlDataFolder='E:\AirwayTreeData';
exportSurfaceFolder='E:\FDsurfaceDataCombo\'; % subject export directory
exportData=0; % export subject datasets
forceupdateXML=0; % force new XML data file creation
forceupdateSurface=0; % force new sufrace extaction from raw XML data

% Set subject and visit variables, settings and selection ranges
load CCStructureSets LUTCCallIDs
subjectIDs=LUTCCallIDs;
subjectIDrange=1; %:length(subjectIDs); 
visitIds={'V1I','V1E','V3I','V3E'};
visitIdrange=1; 

% xml data extraction variables and settings
airwayTreeXMLinputName="ZUNU_vida-xmlTree.xml";
airwayTreeXMLoutputName="AirwayTreeXML.mat";
surfaceDataoutputName="surfaceData_V14.mat";

% segment data extraction settings
segName={'Trachea'}; % ,'LMB','RMB', 'all';
segSurfaces={'inner','outer'}; % Surfaces to process
surfRange=1;

% trachea only settings
segAlignment={'fixed'}; % 'topbot','innerFirst','outerFirst', 'middleThird', 'fixed'

% for rotation purposes use these cutoffs
useTopLung=1; % 1=Top slice selected using top of Lung instead of first completely embedded slice (all perimeter points z>0)
useAroticArch=1; % 1=Bottom slice selected using top of the aortic arch instead of Carina center point (rmb|lmb branch point)
useAlignment=1; % start with sections used for segment alignment
% Segement Slice Selection Range
distanceAbove=99; % distance to use from starting point in mm (+=above -=below) 99=top most embedded slice
distanceBelow=-99; % distance to use from ending point in mm (+=above -=below) -99=Carina


% 2D surface conversion settings
stepSize=0.1; % surface extrapolation X,Y grid spacing in [mm]
shpAlpha=2; % alpha shape value for determining surface perimeter for projection
usenan=1; % preserve nan surface values instead of converting to 0;

% Fractal Variables
maxRange=64; % absolute total deviation allowed in [mm] range=-maxRange/2:maxRange/2
heightAccuracy=1; % final surface Deviation granularlity in [mm]
maxBoxSize=maxRange/heightAccuracy; % total unique grey level variations allowed
sigDecimals=4; 

% ratio ranges to target to esimate for maximum box coverage for values
% where ln(ratio) is approximately split into 0.5 steps 
% value selected is the closest value in each range to the ratio target 
% that maximizes total box coverage of the entire surface. This ensures 
% that ratio values are somewhat evenly distributed for fractal analysis

ratioRanges={2,3,4:5,6:9,10:15,16:25,26:42,43:70,71:115};
ratioTargetL=[2,3,4,7,12,20,33,54,90]; 
ratioTargetU=[2,3,5,8,13,21,34,55,91];

if useTopLung==1||useAroticArch==1
    load('AnatomicalSliceMarkers.mat')
end


%% begin subject data FD extraction

for kSub=subjectIDrange
    for kVis=visitIdrange

        updateXML=forceupdateXML;
        updateSurface=forceupdateSurface;
        segerror=0;

        subjectVisitDataFolder=struct2table(dir(fullfile(importRootFolder,string(subjectIDs(kSub)),visitIds{kVis},'*','*')));

        if size(subjectVisitDataFolder,1)>0 % if subject visit data exists

            disp(strcat(string(subjectIDs(kSub))," Processing visit ",...
                visitIds{kVis}));

            outputSubjectXMLFolder=fullfile(exportXmlDataFolder,string(subjectIDs(kSub)),...
                visitIds{kVis});

            checkFolder(outputSubjectXMLFolder)

            %% Part 1: Raw airway tree xml data extraction
            % checks to see if exising extracted xml data exists and 
            % if not create one and save it for future use

            xmlDataIn=struct2table(dir(fullfile(subjectVisitDataFolder.folder{1},airwayTreeXMLinputName)));
            xmlDataOut=fullfile(outputSubjectXMLFolder,airwayTreeXMLoutputName);

            if size(xmlDataIn,1)==1 && (exist(xmlDataOut,'file')==0 || updateXML==1)

                tic

                airwayTreeXML=readairwayXML(fullfile(xmlDataIn.folder(1,:),xmlDataIn.name(1,:)));
                updateXML=1;

                toc

            elseif exist(xmlDataOut,'file')==2

                xmlfileinfo = matfile(xmlDataOut);

                disp(strcat(string(subjectIDs(kSub))," xml extracted data file for visit ",...
                    visitIds{kVis}," already exists using it instead"));

                % check to see which version of the xmltree data file is 
                % being used and if required update it to the new format.

                if length(who(xmlfileinfo))==1
                    load(xmlDataOut);
                elseif length(who(xmlfileinfo))==6
                    airwayTreeXML=open(xmlDataOut);
                    updateXML=1;
                end

            end

            if exportData==1 && updateXML==1
                save(xmlDataOut,'airwayTreeXML')
            end

            %% Part 2: Segment Point cloud Extraction

            % segment selection
            subjectBaseTable=airwayTreeXML.SegmentNameData(ismember(airwayTreeXML.SegmentNameData.anatomicalName,segName),1:2);

            segIdsUsed=subjectBaseTable.linkIds;
            segIdsList=strcat('SegId_',string(segIdsUsed));

            subjectBaseTable.segId(:,1)=segIdsList;
            subjectBaseTable.id=segIdsUsed;
            subjectBaseTable = movevars(subjectBaseTable, "segId", "Before", "anatomicalName");

            surfaceDataTable=subjectBaseTable;

            segmentNames=fieldnames(airwayTreeXML.SegmentData);
            segmentNames(ismember(string(segmentNames),segIdsList))=[];
            segmentData=rmfield(airwayTreeXML.SegmentData,string(segmentNames));

            updateSurface=1;

            % extract surfaces

            if sum(strcmp(segName,'Trachea'))==1 && strcmp(segAlignment,'fixed')

                tracheaGood=1;

                if strcmp(visitIds{kVis},'V1I')
                    tracheaLocations=V1I_anatomicalSliceLocations(V1I_anatomicalSliceLocations.subjectId==subjectIDs(kSub),:);
                elseif strcmp(visitIds{kVis},'V1E')
                    tracheaLocations=V1E_anatomicalSliceLocations(V1E_anatomicalSliceLocations.subjectId==subjectIDs(kSub),:);
                else
                    % missing key anatomical marker mask
                    tracheaGood=0;
                end

                
                if useTopLung==1 && useAroticArch==1 && tracheaGood==1
                    locTracheaStartSlice=tracheaLocations.topLungSlice;
                    locTracheaEndSlice=tracheaLocations.topArchSlice;
                elseif useTopLung==1 && useAroticArch==0 && tracheaGood==1
                    locTracheaStartSlice=tracheaLocations.topLungSlice;
                    locTracheaEndSlice=0;
                elseif useTopLung==0 && useAroticArch==1 && tracheaGood==1
                    locTracheaStartSlice=0;
                    locTracheaEndSlice=tracheaLocations.topArchSlice;
                elseif (useTopLung==0 && useAroticArch==0) || tracheaGood==0
                    locTracheaStartSlice=nan;
                    locTracheaEndSlice=nan;
                end

            else
                locTracheaStartSlice=nan;
                locTracheaEndSlice=nan;
            end

            if strcmp(segAlignment,'fixed') && (isnan(locTracheaStartSlice) || isnan(locTracheaEndSlice))
                tracheaGood=0;
            end

            tic

            disp('-- Extracting surface data')

            if ~strcmp(segAlignment,'fixed')

                [surfaceMeshs] = extractSurfaces(...
                    segmentData,subjectBaseTable,...
                    airwayTreeXML.voxelDimensions,...
                    'alignment',segAlignment);

            elseif tracheaGood==1 && strcmp(segAlignment,'fixed')

                [surfaceMeshs] = extractSurfaces(...
                    segmentData,subjectBaseTable,...
                    airwayTreeXML.voxelDimensions,...
                    'alignment',segAlignment,...
                    [locTracheaStartSlice,locTracheaEndSlice]);

            else

                disp('unable to process trachea segment using fixed points due to missing value')
                segerror=1;

            end

            toc

            %% Part 3: Generate 2D surfaces from surface meshes for FD analysis
            % check for existing data (to be done later)

            if segerror~=1

                disp('-- Extracting fractal data')


                tic

                surfaceAnalysis.inner=subjectBaseTable;
                surfaceAnalysis.outer=subjectBaseTable;

                for kSeg=1:size(subjectBaseTable,1)

                    if useAlignment==1 % use alignment points trachea mapping cutoff
                        if distanceAbove==99
                            sliceStart=surfaceMeshs.subjectDataTable.outerNonZeroSlice(kSeg);
                        else
                            sliceStart=surfaceMeshs.subjectDataTable.segAlignStart(kSeg);
                        end

                        if distanceBelow==-99
                            sliceEnd=surfaceMeshs.subjectDataTable.sliceCount(kSeg);
                        else
                            sliceEnd=surfaceMeshs.subjectDataTable.segAlignEnd(kSeg);
                        end 
                    else
                        centerLinePointsRawZAxis=surfaceMeshs.(subjectBaseTable.segId(kSeg)).centerLinePointsOriginal(:,3)*airwayTreeXML.voxelDimensions.z;

                    end

                    for kSrf=surfRange

                        % extract segment surfaces

                        ShapeMesh=surfaceMeshs.(segIdsList(kSeg)).(segSurfaces{kSrf}).ShapeMesh;
                        CurveMesh=surfaceMeshs.(segIdsList(kSeg)).(segSurfaces{kSrf}).CurveMesh;
                        TotalMesh=surfaceMeshs.(segIdsList(kSeg)).(segSurfaces{kSrf}).TotalMesh;

                        % trim by slice or by location mm

                        % seperate surface mesh into X,Y (plane), V (distortion) vectors
                        XS=ShapeMesh(:,sliceStart:sliceEnd,1); % same for all surfaces
                        YS=ShapeMesh(:,sliceStart:sliceEnd,2);
                        VS=ShapeMesh(:,sliceStart:sliceEnd,3);

                        XC=CurveMesh(:,sliceStart:sliceEnd,1);
                        YC=CurveMesh(:,sliceStart:sliceEnd,2);
                        VC=CurveMesh(:,sliceStart:sliceEnd,3);

                        XT=TotalMesh(:,sliceStart:sliceEnd,1);
                        YT=TotalMesh(:,sliceStart:sliceEnd,2);
                        VT=TotalMesh(:,sliceStart:sliceEnd,3);

                        % remap surface data onto evenly space grid for fractal
                        % analysis

                        [XqC,YqP,VqC] = surf2grid(XC,YC,VC,stepSize,shpAlpha);
                        [XqT,YqT,VqT] = surf2grid(XT,YT,VT,stepSize,shpAlpha);

                        [XqS,YqS,VqS] = surf2grid(XS,YS,VS,stepSize,shpAlpha);
                        ratiosShape=irdbcRatioSelect(size(VqS,1),size(VqS,2),ratioRanges,ratioTargetL,ratioTargetU);

                        ratiosCurveMap=irdbcRatioSelect(size(VqC,1),size(VqC,2),ratioRanges,ratioTargetL,ratioTargetU);
                        ratiosTotalMap=irdbcRatioSelect(size(VqT,1),size(VqT,2),ratioRanges,ratioTargetL,ratioTargetU);

                        shapeMap=round(VqS./heightAccuracy).*heightAccuracy;
                        shapeZeroth=abs(min(shapeMap(:)))+heightAccuracy;
                        shapeMap=shapeMap+shapeZeroth;

                        curveMap=round(VqC./heightAccuracy).*heightAccuracy;
                        curveZeroth=abs(min(curveMap(:)))+heightAccuracy;
                        curveMap=curveMap+curveZeroth;

                        totalMap=round(VqT./heightAccuracy).*heightAccuracy;
                        totalZeroth=abs(min(totalMap(:)))+heightAccuracy;
                        totalMap=totalMap+totalZeroth;

                        [countShape,ratioShape]=irdbc(shapeMap,ratiosShape,maxBoxSize,shapeZeroth);
                        [countCurve,ratioCurve]=irdbc(curveMap,ratiosCurveMap,maxBoxSize,curveZeroth);
                        [countTotal,ratioTotal]=irdbc(totalMap,ratiosTotalMap,maxBoxSize,totalZeroth);

                        % get linear slope from log(boxRatios)/log(boxCounts) plot of the
                        FD_shape=round(polyfit(log(ratioShape),log(countShape),1),sigDecimals);
                        FD_curve=round(polyfit(log(ratioCurve),log(countCurve),1),sigDecimals);
                        FD_total=round(polyfit(log(ratioTotal),log(countTotal),1),sigDecimals);

                        % get mean square deviation of slope
                        msd_shape=round(sum((log(countShape)-log(ratioShape).*FD_shape(1)-FD_shape(2)).^2)/size(countShape,2),sigDecimals);
                        msd_curve=round(sum((log(countCurve)-log(ratioCurve).*FD_curve(1)-FD_curve(2)).^2)/size(countCurve,2),sigDecimals);
                        msd_total=round(sum((log(countTotal)-log(ratioTotal).*FD_total(1)-FD_total(2)).^2)/size(countTotal,2),sigDecimals);

                        % calculate Surface Roughness
                        SR_shape=round((FD_shape(1)-2)/((FD_shape(1)+2)/2)*2.5*100,sigDecimals);
                        SR_curve=round((FD_curve(1)-2)/((FD_curve(1)+2)/2)*2.5*100,sigDecimals);
                        SR_total=round((FD_total(1)-2)/((FD_total(1)+2)/2)*2.5*100,sigDecimals);

                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).shapeMap=shapeMap; % create surface map w\ shape distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).shapeZeroth=shapeZeroth; % create surface map w\ shape distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).curveMap=curveMap; % create surface map w\ curvature distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).curveZeroth=curveZeroth; % create surface map w\ curvature distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).totalMap=totalMap; % create surface map w\ shape+curvature distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).totalZeroth=totalZeroth; % create surface map w\ shape+curvature distortions

                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).shapeMapRaw=VqS; % create surface map w\ shape distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).curveMapRaw=VqC; % create surface map w\ curvature distortions
                        surfaceMaps.(segIdsList(kSeg)).(segSurfaces{kSrf}).totalMapRaw=VqT; % create surface map w\ shape+curvature distortions

                        % surfaceAnalysis.(segSurfaces{kSrf}).SR_S(kSeg)=SR_shape;
                        % surfaceAnalysis.(segSurfaces{kSrf}).SR_C(kSeg)=SR_curve;
                        % surfaceAnalysis.(segSurfaces{kSrf}).SR_T(kSeg)=SR_total;
                        surfaceAnalysis.(segSurfaces{kSrf}).sliceStart(kSeg)=sliceStart;
                        surfaceAnalysis.(segSurfaces{kSrf}).sliceEnd(kSeg)=sliceEnd;
                        surfaceAnalysis.(segSurfaces{kSrf}).FD_S(kSeg)=FD_shape(1);
                        surfaceAnalysis.(segSurfaces{kSrf}).msd_S(kSeg)=msd_shape;
                        surfaceAnalysis.(segSurfaces{kSrf}).FD_C(kSeg)=FD_curve(1);
                        surfaceAnalysis.(segSurfaces{kSrf}).msd_C(kSeg)=msd_curve;
                        surfaceAnalysis.(segSurfaces{kSrf}).FD_T(kSeg)=FD_total(1);
                        surfaceAnalysis.(segSurfaces{kSrf}).msd_T(kSeg)=msd_total;                      
                        surfaceAnalysis.(segSurfaces{kSrf}).boxRatiosS{kSeg}=ratioShape;
                        surfaceAnalysis.(segSurfaces{kSrf}).boxCountsS{kSeg}=countShape;
                        surfaceAnalysis.(segSurfaces{kSrf}).boxRatiosC{kSeg}=ratioCurve;
                        surfaceAnalysis.(segSurfaces{kSrf}).boxCountsC{kSeg}=countCurve;
                        surfaceAnalysis.(segSurfaces{kSrf}).boxRatiosT{kSeg}=ratioTotal;
                        surfaceAnalysis.(segSurfaces{kSrf}).boxCountsT{kSeg}=countTotal;

                        surfaceAnalysis.(segSurfaces{kSrf}).slopeIntS{kSeg}=FD_shape;
                        surfaceAnalysis.(segSurfaces{kSrf}).slopeIntC{kSeg}=FD_curve;
                        surfaceAnalysis.(segSurfaces{kSrf}).slopeIntSC{kSeg}=FD_total;
                        surfaceAnalysis.(segSurfaces{kSrf}).minDeviation_S(kSeg)=round(min(VqS(:)),2);
                        surfaceAnalysis.(segSurfaces{kSrf}).maxDeviation_S(kSeg)=round(max(VqS(:)),2);
                        surfaceAnalysis.(segSurfaces{kSrf}).minDeviation_C(kSeg)=round(min(VqC(:)),2);
                        surfaceAnalysis.(segSurfaces{kSrf}).maxDeviation_C(kSeg)=round(max(VqC(:)),2);
                        surfaceAnalysis.(segSurfaces{kSrf}).minDeviation_T(kSeg)=round(min(VqT(:)),2);
                        surfaceAnalysis.(segSurfaces{kSrf}).maxDeviation_T(kSeg)=round(max(VqT(:)),2);

                    end


                end
                surfaceAnalysis.surfaceMeshs=surfaceMeshs;
                surfaceAnalysis.surfaceMaps=surfaceMaps;
                surfaceAnalysis.segmentData=segmentData;

                toc

                %% Write out new and updated data files

                if exportData==1 && updateSurface==1

                    outputSubjectSurfaceFolder=fullfile(exportSurfaceFolder,string(subjectIDs(kSub)),...
                    visitIds{kVis});

                    checkFolder(outputSubjectSurfaceFolder)

                    surfaceDataOut=fullfile(outputSubjectSurfaceFolder,surfaceDataoutputName);

                    save(surfaceDataOut,'surfaceAnalysis')

                end

            end
        end

    end
end
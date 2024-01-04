function [airwayTreeXML] = readairwayXML(XMLFileLocation)
% Convert VIDA airway tree xml files into matlab stuctures for processing

        xmlIn=xml2struct(XMLFileLocation);
                
        AirwayBranchPoints=xmlIn.TreeFile.Branchpoints.BP;
        AirwaySegments=xmlIn.TreeFile.Links.Link;
        
        if isfield(xmlIn.TreeFile.SegmentNames,'SegmentName')
            AirwaySegmentsNames=xmlIn.TreeFile.SegmentNames.SegmentName;
        else
            AirwaySegmentsNames=[];
        end
        
        AirwayImageSize=xmlIn.TreeFile.ImageSize.Attributes;
        AirwayVoxelSize=xmlIn.TreeFile.VoxelDimensions.Attributes;
        
        clear xmlIn
        
        %% Segment Name data
        
        SegmentNameData=table;
        
        % build table from character data provide in xml structure
        if ~isempty(AirwaySegmentsNames)
            if length(AirwaySegmentsNames)~=1
                for kName=1:length(AirwaySegmentsNames)
                    
                    if kName==1
                        SegmentNameData=struct2table(structfun(@string,...
                            AirwaySegmentsNames{1, kName}.Attributes,...
                            'UniformOutput',false));

                    else
                        nextRow=struct2table(structfun(@string,...
                            AirwaySegmentsNames{1, kName}.Attributes,...
                            'UniformOutput',false));
                        SegmentNameData=[SegmentNameData;nextRow];
                    
                    end
                end
                
            elseif length(AirwaySegmentsNames)==1
                
                SegmentNameData=struct2table(structfun(@string,...
                    AirwaySegmentsNames.Attributes,...
                    'UniformOutput',false));
                
            end
            
            SegmentNameData=movevars(SegmentNameData,...
                'linkIds','before',2);
            SegmentNameData=movevars(SegmentNameData,...
                'startBpId','before',3);
            
            % convert imported table strings into numerial doubles
            for kCol=2:size(SegmentNameData,2)
                SegmentNameData.(kCol)=str2double(SegmentNameData{:,kCol});
            end
        end
        % done branch point data extraction
        
        % branchPoint Voxel Locations
        BranchPointXYZData=table;

        % build table from character data provide in xml structure
        for kBranch=1:length(AirwayBranchPoints)
            if kBranch==1
                
                BranchPointXYZData=struct2table(structfun(@string,...
                    AirwayBranchPoints{1, kBranch}.Attributes,...
                    'UniformOutput',false));
            else
                nextRow=struct2table(structfun(@string,...
                    AirwayBranchPoints{1, kBranch}.Attributes,...
                    'UniformOutput',false));

                BranchPointXYZData=[BranchPointXYZData;nextRow];
            end
        end
        
        % convert imported table strings into numerial doubles
        for kCol=1:size(BranchPointXYZData,2)

            BranchPointXYZData.(kCol)=...
                str2double(BranchPointXYZData{:,kCol});
        end
        % done branch point data extraction
        
        % Segment and connecting point Voxel data
        BranchSegmentData=table;
        SegmentData=struct;
        % ConnectedPointXYZ=struct;
        
        % build table from character data provide in xml structure
        for kSegment=1:length(AirwaySegments)
            
            ConnectedPointData=table;
            
            if kSegment~=1
                
                nextRow=struct2table(structfun(@string,...
                    AirwaySegments{1, kSegment}.Attributes,...
                    'UniformOutput',false));
                
                % check for missing variable columns Master Table
                MissingVar=setdiff(nextRow.Properties.VariableNames,...
                    BranchSegmentData.Properties.VariableNames);
                
                BranchSegmentData = [BranchSegmentData,...
                    array2table(nan(height(BranchSegmentData),...
                    numel(MissingVar)), 'VariableNames', MissingVar)];
                
                % check for missing variable columns new row
                MissingVar=setdiff(...
                    BranchSegmentData.Properties.VariableNames, ...
                    nextRow.Properties.VariableNames);
                
                nextRow = [nextRow, ...
                    array2table(nan(height(nextRow), ...
                    numel(MissingVar)), 'VariableNames', MissingVar)];
                
                % merge next row with previous table
                BranchSegmentData=[BranchSegmentData;nextRow];
                
            else
                
                BranchSegmentData=struct2table(structfun(@string,...
                    AirwaySegments{1, kSegment}.Attributes,...
                    'UniformOutput',false));
                BranchSegmentData=movevars(BranchSegmentData,...
                    'id','before',1);
                
            end
            
            % Connecting Points
            % make sure CP field exists
            if isfield(AirwaySegments{1, kSegment},'CP')==1 
                
                for kConnectedPoints=1:length(...
                        AirwaySegments{1, kSegment}.CP)
                    
                    % Build ConnectedPointData Table
                    
                    if kConnectedPoints~=1 && ...
                            isstruct(AirwaySegments{1, kSegment}.CP)==0
                        
                      nextRow=struct2table(structfun(@string,...
                      AirwaySegments{1, kSegment}.CP{1,...
                      kConnectedPoints}.Attributes,...
                            'UniformOutput',false));
                        
                        % check for missing variable columns Master Table
                        MissingVar=setdiff(...
                            nextRow.Properties.VariableNames,...
                            ConnectedPointData.Properties.VariableNames);
                        ConnectedPointData = [ConnectedPointData ...
                            array2table(nan(height(ConnectedPointData), ...
                            numel(MissingVar)), 'VariableNames',...
                            MissingVar)];
                        
                        % check for missing variable columns new row
                        MissingVar=setdiff(...
                            ConnectedPointData.Properties.VariableNames,...
                            nextRow.Properties.VariableNames);
                        nextRow = [nextRow ...
                            array2table(nan(height(nextRow), ...
                            numel(MissingVar)), 'VariableNames',...
                            MissingVar)];
                        
                        % merge next row with previous table
                        ConnectedPointData=[ConnectedPointData;nextRow];
                        
                    elseif isstruct(AirwaySegments{1, kSegment}.CP)==0
                        
                        ConnectedPointData=...
                            struct2table(structfun(@string,...
                            AirwaySegments{1, kSegment}.CP{1,...
                            kConnectedPoints}.Attributes,...
                            'UniformOutput',false));
                        ConnectedPointData=movevars(ConnectedPointData,...
                            'id','before',1);
                    end
                    
                    % if Connecting points exist
                    if (isstruct(AirwaySegments{1, kSegment}.CP)==0 && ... 
                            isfield(AirwaySegments{1, kSegment}.CP{1,...
                            kConnectedPoints},'Text')==0) 
                        %if no Text Variable
                        %Inner and Outerborder points exist
                        
                        % connectedPointsValues
                        
                        innerBorderPoints=...
                            str2num(AirwaySegments{1, kSegment}.CP{1,...
                            kConnectedPoints}.InnerBorderPnts.Text);
                        innerBorderPoints=reshape(innerBorderPoints,...
                            [3,length(innerBorderPoints)/3])';
                        innerBorderPoints=array2table(innerBorderPoints,...
                            'VariableNames',{'x','y','z'});
                        
                        outerBorderPoints=str2num(...
                            AirwaySegments{1, kSegment}.CP{1,...
                            kConnectedPoints}.OuterBorderPnts.Text);
                        outerBorderPoints=reshape(outerBorderPoints,...
                            [3,length(outerBorderPoints)/3])';
                        outerBorderPoints=array2table(outerBorderPoints,...
                            'VariableNames',{'x','y','z'});
                        
                        SegmentData.(strcat('SegId_',...
                            BranchSegmentData.id(kSegment))). ...
                            ConnectedPointXYZ.(strcat('CPtId_',...
                            ConnectedPointData.id(kConnectedPoints))).inner...
                            =innerBorderPoints;
                        SegmentData.(strcat('SegId_',...
                            BranchSegmentData.id(kSegment))). ...
                            ConnectedPointXYZ.(strcat('CPtId_',...
                            ConnectedPointData.id(kConnectedPoints))).outer...
                            =outerBorderPoints;
                    end
                end
                
                for kCol=1:size(ConnectedPointData,2)
                    ConnectedPointData.(kCol)=...
                        str2double(ConnectedPointData{:,kCol});
                end
                
                SegmentData.(strcat('SegId_',...
                    BranchSegmentData.id(kSegment))). ...
                    ConnectedPointData=ConnectedPointData;
            end
            
            
        end
        
        % convert imported numerical string values into numerical doubles
        
        strCols=["direction","subLobe","lobe"];
        
        for kCol=1:size(BranchSegmentData,2)
            if ~ismember(...
                    BranchSegmentData.Properties.VariableNames{kCol},...
                    strCols)
                BranchSegmentData.(kCol)=...
                    str2double(BranchSegmentData{:,kCol});
            end
        end
        
        % Lung Slice size and Voxel Dimensions
        
        imageDimensions=struct2table(AirwayImageSize);
        voxelDimensions=struct2table(AirwayVoxelSize);
        
        for kCol=1:3
            imageDimensions.(kCol)=str2double(imageDimensions{:,kCol});
            voxelDimensions.(kCol)=str2double(voxelDimensions{:,kCol});
        end
              
        airwayTreeXML.imageDimensions=imageDimensions;
        airwayTreeXML.voxelDimensions=voxelDimensions;
        airwayTreeXML.BranchPointXYZData=BranchPointXYZData;
        airwayTreeXML.BranchSegmentData=BranchSegmentData;
        airwayTreeXML.SegmentData=SegmentData;
        airwayTreeXML.SegmentNameData=SegmentNameData;

end
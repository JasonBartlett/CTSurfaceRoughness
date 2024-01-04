function [ratios2use] = irdbcRatioSelect(M,N,ratioRanges,ratioTargetL,ratioTargetU)
%irdbcRatioSelect Selects ratios to use based on maximizing box coverage
%   selects ratios based log value ranges such that one value is selected
%   from each ratio range that provides maximum box coverage

ratios2use=[];

for kRat=1:length(ratioRanges)

    ratios2consider=ratioRanges{kRat};
    boxSizesM=floor(M./ratios2consider);
    boxSizesN=floor(N./ratios2consider);

    boxCoverageM=boxSizesM.*ratios2consider-M;
    boxCoverageN=boxSizesN.*ratios2consider-N;

    boxCoverageSumMN=boxCoverageM+boxCoverageN;

    ratiosSelected=ratios2consider(abs(boxCoverageSumMN)==(min(abs(boxCoverageSumMN))));

    while length(ratiosSelected)>1

        % if there are multiple minimal coverage
        % candidates pick the once closest to expected
        % target value

        ratiosSelected=ratiosSelected(abs(ratiosSelected-ratioTargetL(kRat))...
            ==min(abs(ratiosSelected-ratioTargetL(kRat))));

        % if both values are equal distance away then
        % select using target 2

        ratiosSelected=ratiosSelected(abs(ratiosSelected-ratioTargetU(kRat))...
            ==min(abs(ratiosSelected-ratioTargetU(kRat))));
    end

    ratios2use=unique([ratios2use,ratiosSelected]);

end

end
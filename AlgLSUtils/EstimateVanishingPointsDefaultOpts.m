function opts = EstimateVanishingPointsDefaultOpts(method)

% generate default options
% $Id: EstimateVanishingPointsDefaultOpts.m 1262 2012-01-20 01:55:00Z faraz $

% options
opts.conditioning = false;
opts.outputFormat = 2;

opts.kActionFunctionIndex = 2;
opts.kNormalizationTol = 1e-16;
opts.kPolyResidualTol = inf;
opts.kSolutionInClassTol = .99;
opts.kMulMatrixCondToAbort = 1e14;
switch method
    case 'minimal'
        opts.Ne = 4;
        opts.NeMax = 4;
    case 'relaxed'
        opts.Ne = 7;
        opts.NeMax = 7;        
    case 'nonrelaxed'
        opts.Ne = 14;
        opts.NeMax = 14;        
    otherwise
        error('Uknown method requested for estimating vanishing points.');
end

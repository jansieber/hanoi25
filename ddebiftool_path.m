function ddebiftool_path(base)
if nargin<1
    base=[pwd(),filesep,'DDE_Biftool2025',filesep];
end
addpath([base,filesep,'ddebiftool'],...
    [base,filesep,'ddebiftool_extra_psol'],...
    [base,filesep,'ddebiftool_utilities'],...
    [base,filesep,'ddebiftool_extra_rotsym'],...
    [base,filesep,'ddebiftool_extra_nmfm'],...
    [base,filesep,'ddebiftool_extra_symbolic'],...
    [base,filesep,'ddebiftool_coco']);
end
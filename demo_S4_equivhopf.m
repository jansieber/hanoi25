%%
clear;
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
format short g
load('S4_demo_psol_results.mat');
%% Find a few symmetries that one can branch off into
% The function check_hopf below should return 2 zero singular value to
% satisfy the equivariant Hopf bifurcation theorem
check_hopf=@(cnd)SetupHopf(funcs,eq_1234ref,eq_1234bifind,...
    'initcond',cnd,'usercond',cnd,'output','sv','outputfuncs',true);
%% Permutation (1234)_0 (time shift 0)
ev1234_0=permfix('stst','p1234_1_4',n_osc,'v','perm',[1,2,3,4],'rotation',[0,1]);
sv_1234_0=check_hopf(ev1234_0)'
%% Permutation (12)_0 only is not enough: we have four zero singular values for time shifts
%0 or 1/2. Time shift 1/4 would be ok
ev12=permfix('stst','p12',n_osc,'v','perm',[1,2],'rotation',[0,1]);
sv_12=check_hopf(ev12)'
%% Permutation (1234)_{1/4} (time shift 1/4)
ev1234_1_4=permfix('stst','p1234_1_4',n_osc,'v','perm',[1,2,3,4],'rotation',[1,4]);
sv_1234_1_4=check_hopf(ev1234_1_4)'
%% Permutation  (123)_0 (time shift 0)
ev123=permfix('stst','pv123',n_osc,'v','perm',[1,2,3],'rotation',[0,1]);
sv_123=check_hopf(ev123)'
%% Permutation (12)_0(34)_0:  1=2 (time shift 0),3=4 (time shift 0), 1=3 (time shift 1/2)
ev12_34=[permfix('stst','pv12_0_1',n_osc,'v','perm',[1,2],'rotation',[0,1]),...
         permfix('stst','pv34_0_1',n_osc,'v','perm',[3,4],'rotation',[0,1])];
sv_12_34=check_hopf(ev12_34)'
%% Permutation (12)_{1/2}(34)_0: 1=2 (time shift 1/2),3=4 (time shift 0)
ev12s_34=[permfix('stst','p12_1_2',n_osc,'v','perm',[1,2],'rotation',[1,2]),...
          permfix('stst','p34_0_1',n_osc,'v','perm',[3,4],'rotation',[0,1])];
sv_12s_34=check_hopf(ev12s_34)';
%% Permutation (12)_{1/4}: 1=2 (time shift 1/4) is sufficient
ev12_1_4=permfix('stst','p12',n_osc,'v','perm',[1,2],'rotation',[1,4]);
sv_12_1_4=check_hopf(ev12_1_4)'

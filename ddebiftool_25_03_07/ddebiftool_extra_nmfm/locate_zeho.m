function [out,success,zhfuncs]=locate_zeho(funcs,branch,ind,varargin)
if funcs.tp_del
    sys_tau=@(it,x,p)sys_tau_zeho(it,x,p,funcs);
else
    sys_tau=funcs.sys_tau;
end
zhcond=dde_sys_cond_create('fun',@(p,pref)sys_cond_zeho(p,pref,@get_comp),'name','sys_cond_zeho','reference',true);
zhfuncs=set_funcs(...
    'sys_rhs',@(x,p)sys_rhs_zeho(funcs,0,x,p),...
    'sys_tau',sys_tau,...
    'sys_dirderi',{@(x,p,dx,dp)sys_rhs_zeho(funcs,1,x,p,dx,dp)},...
    'sys_ntau',funcs.sys_ntau,'sys_cond',zhcond,...
    'x_vectorized',true,'p_vectorized',true);
zhfuncs.get_comp=@get_comp;
zhfuncs.userfuncs=funcs;
branch=replace_branch_pars(branch,branch.parameter.free,varargin);
zehoini=branch.point(ind);
zhmth=setfield(branch.method.point,'postprocess','');
zhmth.extra_condition=true;
zeho0=zehoext_from_stst(funcs,zehoini,branch.method.stability);
[zeho,success]=p_correc(zhfuncs,zeho0,...
    [branch.parameter.free,length(zeho0.parameter)],[],zhmth,[],zeho0);
if ~success
    warning('locate_zeho:notconv','locate_zeho: not converged');
end
out=get_comp(zeho,'solution');
switch zehoini.kind
    case 'fold'
        out=dde_fold_create('point',out,'v',get_comp(zeho,'q0'),'flag','zeho');
    case 'hopf'
        out=dde_hopf_create('point',out,'v',get_comp(zeho,'q1'),...
            'omega',get_comp(zeho,'omega'),'flag','zeho');
end
out.nvec.q0=get_comp(zeho,'q0');
out.nvec.q1=get_comp(zeho,'q1');
out.nvec.omega=get_comp(zeho,'omega');
out.stability=p_stabil(funcs,out,branch.method.stability);
end
%%
function out=get_comp(zeho,comp)
        n=length(zeho.x)/4;
switch comp
    case 'solution'
        out=dde_stst_create('x',real(zeho.x(1:n)),'parameter',real(zeho.parameter(1:end-1)));
    case 'x'
        out=zeho.x(1:n);
    case 'ind_x'
        out=1:n;
    case 'parameter'
        out=zeho.parameter(1:end-1);
    case 'ind_parameter'
        out=1:length(zeho.parameter)-1;
    case 'omega'
        out=zeho.parameter(end);
    case 'q0'
        out=zeho.x(n+(1:n));
    case 'ind_q0'
        out=n+(1:n);
    case 'q1'
        out=zeho.x(2*n+(1:n))+1i*zeho.x(3*n+(1:n));
    case 'ind_rq1'
        out=2*n+(1:n);
    case 'ind_iq1'
        out=3*n+(1:n);
end
end
%%
function [zehoext,zeho0]=zehoext_from_stst(funcs,stst,stmth)
if ~isfield(stst,'stability')|| isempty(stst.stability)
    stst.stability=p_stabil(funcs,stst,stmth);
end
if strcmp(stst.kind,'stst')
    zeho=p_tofold(funcs,stst);
else
    zeho=stst;
end
zeho=rmfield(zeho,'nvec');
zeho0=nmfm_zeho(funcs,p_tozeho(zeho));
zehoext=dde_stst_create('x',[zeho0.x;zeho0.nvec.q0;real(zeho0.nvec.q1);imag(zeho0.nvec.q1)],...
    'parameter',[zeho0.parameter,zeho0.nvec.omega]);
end
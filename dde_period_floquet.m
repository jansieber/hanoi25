function stability=dde_period_floquet(Afun,delays,period,varargin)
%% compute stability information for stst 
% INPUT:
%   funcs problem functions
%	pt solution point
%	method method parameters 
% OUTPUT:
%	stability stability information
%
default={'check',true,'geteigenfuncs',false,'fill',[],'closest',[],...
    'eigmatrix','full','delay_accuracy',1e-8,'collocation_parameters',[],...
    'max_number_of_eigenvalues',20,'minimal_modulus',0,'degree',4,'intervals',20};
[opts,pass_on]=dde_set_options(default,varargin,'pass_on');
nd=length(delays);
bmesh=dde_coll_meshfill(linspace(0,1,opts.intervals+1),opts.degree,'purpose','storage','grid','cheb');
cmesh=dde_coll_meshfill(linspace(0,1,opts.intervals+1),opts.degree,'purpose','collocation','grid','gauss');
[ntb,ntc]=deal(length(bmesh),length(cmesh));
At=Afun(cmesh*period);
nx=size(At,1);
Amat=mshape(At,'r',{nx,nx*nd,length(cmesh)},'apply',@sparse_blkdiag,'state','sparse','struct',...
    'r',{nx,ntc,nx,nd,ntc},'p',[1,2,3,5,4],'r',{nx*ntc,nx*ntc*nd},'mat');
pt=dde_psol_create('profile',zeros(nx,ntb),'mesh',bmesh,'degree',opts.degree,'period',period);
tc=(cmesh(ones(nd,1),:)-repmat(delays(:)/period,1,ntc))';
[~,J0,ptext]=dde_coll_wrap_eva(pt,tc(:).','wrapJ',false,'kron',true);
J1=dde_coll_eva(pt,cmesh,'diff',1,'output','matrix','kron',true);
M=Amat*J0;
em=repmat(ptext.mesh,nx,1);
M(:,em(:)>=0)=M(:,em(:)>=0)-J1/period;
%% compute eigenvalues and eigenfunctions
[Margs,MPfun]=dde_psol_monodromy(M,opts.closest,opts.eigmatrix);
isef=opts.geteigenfuncs;
[s,ef]=feval(['dde_psol_mult_app_',opts.eigmatrix],Margs,...
    'geteigenfuncs',isef,'closest',opts.closest',...
    'max_number_of_eigenvalues',opts.max_number_of_eigenvalues,pass_on{:});
assert(size(ef,2)==length(s)|~isef);
%% remove NaNs and Infs
isval=isfinite(s);
s=s(isval);
ef=ef(:,isval&isef);
%% sort, ensuring that positive imaginary parts come first
[~,im]=sort(imag(s),'descend');
s=s(im);
if isef
    ef=ef(:,im);
end
if isempty(opts.closest)
    [~,is]=sort(abs(s),'descend');
else
    [~,is]=sort(abs(s-opts.closest(1)),'ascend');
end
s=s(is);
if isef
    ef=ef(:,is);
end
%% remove unwanted eigenvalues   
sel= abs(s)>opts.minimal_modulus & (1:length(s))'<=opts.max_number_of_eigenvalues;
s=s(sel);
ef=ef(:,sel&isef);
neig=length(s);
err=NaN(1,neig);
if isef
    %% expand (if necesssary) eigenfunction to period interval + delay time
    ef=MPfun(ef);
    if opts.check==1
        y=reshape(ef,nx,size(em,2),neig);
        exsel=ptext.mesh<=ptext.mesh(end)-1;
        exm1=ptext.mesh(exsel);
        for i=1:neig
            y1=dde_coll_eva(y(:,:,i),ptext.mesh,exm1+1,opts.degree);
            y0=y(:,exsel,i);
            err(i)=max(norm(M*ef(:,i),'inf'),norm(y0(:)*s(i)-y1(:),'inf'));
        end
    else
        err=NaN(1,neig);
    end
    %% convert eigenvectors to point structures & add interpolation error
    ef=ef(em>=0&em<=1,:);
    efp=dde_point_from_x([ef;zeros(1,neig)+pt.period],pt,[]);
    for i=1:neig % scale to max norm =1
        efp(i).profile=efp(i).profile/max(abs(efp(i).profile(:)));
    end
    if opts.check
        for i=1:neig
            err(i)=max(err(i),dde_coll_error(efp(i),'output','norm'));
        end
    end
else
    efp=repmat(pt,1,0);
end
%% fill entries in s with 'fill'
if ~isempty(opts.fill)
    s=cat(1,s,repmat(opts.fill,opts.max_number_of_eigenvalues-length(s),1));
end
%%
stability=struct('mu',s,'eigenfuncs',efp,'err',err);
end

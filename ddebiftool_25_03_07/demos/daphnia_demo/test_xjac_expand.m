[fullargs,dimx,nvecx]=arg_flatten(fmt,arg_exp({1,'I'},0,0));
Jx=dde_coll_rhs(funcs,fullargs{:});
Jxr=reshape(Jx,[size(Jx,1),dimx]);
Jxo=Jc{ia.x};
%%
ir=find(abs(Jxo(:)-Jxr(:))>1e-10);
[jf,jx,jt,jd,jntst]=ind2sub(size(Jxr),ir);
dd=diff([jf,jx,jt,jd,jntst]',[],2);
Jco=cellfun(@(a,f){reshape(a,[size(a,1),f])},mat2cell(Jall,size(Jall,1),[nvec{:}]),dim);
Jxo0=Jco{ia.x};
// Isolated partial-universe BM / monotone-within-live-component prototype.
// The host namespace is Matchbox; all methods belong to FastMatchingEngineV2.
// Reuse an entry's leaf-to-ancestor route while descending into its first
// child. Child-base identities find the base child without another ancestry
// walk. Every other recursive branch has a disjoint child subtree.
inline void FastMatchingEngineV2::supportPath(
    Size blossom,Size entry,std::vector<Size>& out)
{
  struct Task {Size b,x,route,pos;bool reverse;};
  std::vector<std::vector<Size>> routes;
  std::vector<Task> tasks(1,Task{blossom,entry,cNullItm,0,false});out.clear();
  while(!tasks.empty()) {
    Task t=tasks.back();tasks.pop_back();
    if(t.b>=mBlsm.size()||!mBlsm[t.b].alive||t.x>=mNumVtxs)
      throw std::runtime_error("support path invalid task");
    if(trivial(t.b)||t.x==mBlsm[t.b].base) {
      if(trivial(t.b)&&t.b!=t.x)throw std::runtime_error("support path wrong atomic entry");
      out.push_back(t.x);continue;
    }
    if(t.route==cNullItm) {
      t.route=routes.size();routes.push_back(std::vector<Size>());
      Size b=t.x;
      for(Size guard=0;;++guard) {
        if(b>=mBlsm.size()||guard>mBlsm.size())throw std::runtime_error("support entry route escapes blossom");
        routes[t.route].push_back(b);MS_ADD(routeLinks,1);
        if(b==t.b)break;b=mBlsm[b].parent;
      }
      t.pos=routes[t.route].size()-1;
    }
    if(!t.pos||routes[t.route][t.pos]!=t.b)throw std::runtime_error("support cached route mismatch");
    const Size entryChild=routes[t.route][t.pos-1];
    const auto& B=mBlsm[t.b];const Size k=B.children.size();
    Size ib=cNullItm,cur=cNullItm;
    for(Size j=0;j<k;++j) {
      MS_ADD(childInspections,1);
      if(B.children[j]==entryChild)cur=j;
      if(mBlsm[B.children[j]].base==B.base)ib=j;
    }
    if(!k||B.childEdges.size()!=k||ib==cNullItm||cur==cNullItm)
      throw std::runtime_error("support path child/base identity missing");
    std::vector<std::pair<Size,Size>> chain(1,std::make_pair(entryChild,t.x));
    for(Size guard=0;cur!=ib;++guard) {
      if(guard>k)throw std::runtime_error("support child cycle");
      const Size rel=(cur+k-ib)%k;Size by,bx,edge;
      if(!(rel&1)){by=(cur+k-1)%k;bx=(cur+k-2)%k;edge=B.childEdges[bx];cur=bx;}
      else {bx=(cur+1)%k;by=(bx+1)%k;edge=B.childEdges[bx];cur=by;}
      const Size x=edgeSrc(edge),y=edgeDst(edge);
      if(x>=mNumVtxs||y>=mNumVtxs)throw std::runtime_error("support child edge invalid");
      if(!(rel&1)) {
        chain.push_back(std::make_pair(B.children[by],y));chain.push_back(std::make_pair(B.children[bx],x));
      } else {
        chain.push_back(std::make_pair(B.children[bx],x));chain.push_back(std::make_pair(B.children[by],y));
      }
    }
    auto push=[&](Size j){
      tasks.push_back(Task{chain[j].first,chain[j].second,j?cNullItm:t.route,j?0:t.pos-1,(bool)(t.reverse^(bool)(j&1))});
    };
    if(!t.reverse){for(Size j=chain.size();j>0;--j)push(j-1);}
    else {for(Size j=0;j<chain.size();++j)push(j);}
  }
  if(out.empty()||out.front()!=entry||out.back()!=mBlsm[blossom].base||!(out.size()&1))
    throw std::runtime_error("support path boundary/parity mismatch");
  MS_ADD(paths,1);MS_ADD(pathAtoms,out.size());MS_PEAK(maxPath,out.size());
}

inline void FastMatchingEngineV2::supportAdd(
    Size v,Size parent,bool singleton,int kind)
{
  if(v>=mNumVtxs||parent>mNumVtxs||!mInT[parent]||mInT[v]||mBM.blossom_of(v)!=cNullItm)
    throw std::runtime_error("support add is not a new unowned leaf");
  const Size tree=mSupportTreeOf[parent];
  if(tree>=mNumVtxs)throw std::runtime_error("support add lacks live component");
  mAuxNCA.addLeaf(v,parent);
#ifdef MBX_GABOW_SEARCH_AUDIT
  audNoteAddLeaf(v);
#endif
  mInT[v]=1;mBM.set_depth(v,mAuxNCA.depth(v));
  mSupportTreeOf[v]=tree;mSupportNext[v]=mSupportHead[tree];mSupportHead[tree]=v;
  if(singleton){mSupportSingleton[0]=v;mBM.make_blossom(mSupportSingleton);}
  if(kind==0)MS_ADD(innerAdds,1);else if(kind==1)MS_ADD(outerAdds,1);else MS_ADD(shrinkAdds,1);
}

inline Size FastMatchingEngineV2::supportInstallPath(Size b,Size entry,Size anchor)
{
  if(anchor>=mNumVtxs||!mInT[anchor])throw std::runtime_error("inner support anchor absent");
  std::vector<Size> path;supportPath(b,entry,path);Size prev=anchor;
  for(Size v:path) {
    if(!mInT[v])supportAdd(v,prev,true,0);
    else if(mAuxNCA.parent(v)!=prev||mSupportTreeOf[v]!=mSupportTreeOf[anchor])
      throw std::runtime_error("inner support path would reparent a live vertex");
    prev=v;
  }
  return prev;
}

inline void FastMatchingEngineV2::supportInner(Size b)
{
  const Size pred=mBlsm[b].pred;
  if(pred==cNullItm)throw std::runtime_error("inner support has no entering edge");
  supportInstallPath(b,edgeDst(pred),edgeSrc(pred));
}

inline Size FastMatchingEngineV2::supportMergePath(Size from,Size ancestor)
{
  if(from>=mNumVtxs||ancestor>=mNumVtxs||!mInT[from]||!mInT[ancestor]||
     mSupportTreeOf[from]!=mSupportTreeOf[ancestor])
    throw std::runtime_error("support merge crosses live components");
  for(Size guard=0;guard<=mNumVtxs;++guard) {
    const Size owner=mBM.blossom_of(from),target=mBM.blossom_of(ancestor);
    if(owner==cNullItm||target==cNullItm)throw std::runtime_error("support merge missing owner");
    if(owner==target)return owner;
    const Size root=mBM.mBlossoms[owner].rootVertex,parent=mAuxNCA.parent(root);
    if(parent>=mNumVtxs||!mInT[parent]||mAuxNCA.depth(parent)<mAuxNCA.depth(ancestor))
      throw std::runtime_error("support merge misses ancestor");
    const Size up=mBM.blossom_of(parent);
    if(up==cNullItm||up==owner)throw std::runtime_error("support owner root is not a boundary");
    mBM.merge(owner,up,&mMergeDiscarded);++mProfMergeCount;MS_ADD(merges,1);from=parent;
  }
  throw std::runtime_error("support merge cycle");
}

inline void FastMatchingEngineV2::supportRepublish()
{
  GabowDeferredTwins saved;
  for(Size pl:mMergeDiscarded) {
    const Size u=pl/mAdjStride,i=pl%mAdjStride;
    if(u>=mNumVtxs||i>=mAdj[u].size())throw std::runtime_error("support discard invalid payload");
    stageGabowTwins(u,i,saved);
    mD3QueuedPopulation.set(mD3Queued,u,i,(char)0);
    mBmHandleU[u][i]=cNullItm;mBmHandleV[u][i]=cNullItm;
  }
  removeStagedGabowTwins(saved);
  for(Size pl:mMergeDiscarded)pushD3(pl/mAdjStride,pl%mAdjStride);
  mMergeDiscarded.clear();
}

inline void FastMatchingEngineV2::supportOuter(Size b)
{
  std::vector<Size> collected;MembershipView vertices=&mBlsm[b].vertices;
  if(vertices->empty()){collectBlossomVertices(b,&collected);vertices=&collected;}
  const Size base=mBlsm[b].base,pred=mBlsm[b].pred;
  if(base>=mNumVtxs||pred==cNullItm||edgeDst(pred)!=base)
    throw std::runtime_error("outer support matched-boundary mismatch");
  mMergeDiscarded.clear();
  if(!mInT[base]) {
    // Ordinary grow: the entire new outer blossom is absent. A single bulk
    // make_blossom avoids one merge/heap-maintenance operation per new atom.
    for(Size v:*vertices)if(mInT[v]||mBM.blossom_of(v)!=cNullItm)
      throw std::runtime_error("fresh outer support overlaps existing universe");
    supportAdd(base,edgeSrc(pred),false,1);
    for(Size v:*vertices)if(v!=base)supportAdd(v,base,false,1);
    if(vertices->size()==1){mSupportSingleton[0]=base;mBM.make_blossom(mSupportSingleton);}
    else {std::vector<Size> atoms(vertices->begin(),vertices->end());mBM.make_blossom(atoms);}
  } else {
    // Expand: preserve the old path. Add only absent atoms as new leaves,
    // then coalesce the new outer set along actual support-parent links.
    for(Size v:*vertices)if(!mInT[v])supportAdd(v,base,true,1);
    for(Size v:*vertices)supportMergePath(v,base);
  }
  const Size owner=mBM.blossom_of(base);mMatchboxToBm[b]=owner;
  for(Size v:*vertices){
    if(mBM.blossom_of(v)!=owner)throw std::runtime_error("new outer support ownership incomplete");
    mMatchboxToBm[v]=owner;
  }
  supportRepublish();
}

inline void FastMatchingEngineV2::supportBlossom(Size b,Size u,Size v)
{
  if(!mInT[u]||!mInT[v]||mSupportTreeOf[u]!=mSupportTreeOf[v])
    throw std::runtime_error("blossom support endpoints are not in one component");
  mMergeDiscarded.clear();
  for(Size w:mBlsm[b].vertices)if(!mInT[w]) {
    supportAdd(w,u,true,2);
    mBM.merge(mBM.blossom_of(w),mBM.blossom_of(u),&mMergeDiscarded);++mProfMergeCount;MS_ADD(merges,1);
  }
  const Size ancestor=mAuxNCA.nca(u,v);
  if(ancestor>=mNumVtxs)throw std::runtime_error("same-tree blossom support has dummy NCA");
  supportMergePath(u,ancestor);supportMergePath(v,ancestor);
  const Size owner=mBM.blossom_of(ancestor);mMatchboxToBm[b]=owner;
  for(Size child:mBlsm[b].children)mMatchboxToBm[child]=owner;
  for(Size w:mBlsm[b].vertices){
    if(mBM.blossom_of(w)!=owner)throw std::runtime_error("blossom support cycle not fully coalesced");
    mMatchboxToBm[w]=owner;
  }
  supportRepublish();
}

inline void FastMatchingEngineV2::supportRetireTree(Size tree)
{
  if(tree>=mSupportHead.size()||mSupportHead[tree]==cNullItm)return;
  if(++mSupportOwnerEpoch==0){std::fill(mSupportOwnerStamp.begin(),mSupportOwnerStamp.end(),0);mSupportOwnerEpoch=1;}
  std::vector<Size> owners;Size v=mSupportHead[tree];
  while(v!=cNullItm) {
    if(v>=mNumVtxs||mSupportTreeOf[v]!=tree)throw std::runtime_error("support retirement list corrupt");
    const Size next=mSupportNext[v],owner=mBM.blossom_of(v);
    if(owner==cNullItm)throw std::runtime_error("support retirement atom has no owner");
    if(owner>=mSupportOwnerStamp.size())mSupportOwnerStamp.resize(owner+1,0);
    if(mSupportOwnerStamp[owner]!=mSupportOwnerEpoch){mSupportOwnerStamp[owner]=mSupportOwnerEpoch;owners.push_back(owner);}
    retireGabowForestVertex(v);mSupportTreeOf[v]=cNullItm;mSupportNext[v]=cNullItm;
    mMatchboxToBm[v]=cNullItm;MS_ADD(retiredAtoms,1);v=next;
  }
  mSupportHead[tree]=cNullItm;mMergeDiscarded.clear();
  for(Size owner:owners){mBM.dispose_blossom(owner,&mMergeDiscarded);++mProfDisposeCount;MS_ADD(retiredOwners,1);}
  // All consumed owners are fully gone before generation-guarded twin
  // cleanup and refresh. Other trees' forests/owners are not reset.
  supportRepublish();MS_ADD(retiredTrees,1);
}

inline void FastMatchingEngineV2::supportVerifyStable(const char* where) const
{
#ifdef MBX_MONOTONE_SUPPORT_VERIFY
  if(mBackend!=kGabow2018)return;
  std::vector<char> listed(mNumVtxs,0);
  for(Size tree=0;tree<mSupportHead.size();++tree) {
    Size steps=0;
    for(Size v=mSupportHead[tree];v!=cNullItm;v=mSupportNext[v]) {
      if(v>=mNumVtxs||listed[v]||!mInT[v]||mSupportTreeOf[v]!=tree||++steps>mNumVtxs)
        throw std::runtime_error(std::string("support list invariant at ")+where);
      listed[v]=1;
    }
  }
  for(Size v=0;v<mNumVtxs;++v) {
    Size b=v,steps=0;
    while(mBlsm[b].parent!=cNullItm){b=mBlsm[b].parent;if(++steps>mBlsm.size())throw std::runtime_error("support logical parent cycle");}
    const bool inSearch=mBlsm[b].status==BS_EVEN||mBlsm[b].status==BS_ODD;
    if(mBlsm[b].status==BS_EVEN&&!mInT[v])throw std::runtime_error("outer atom missing from support");
    if(!inSearch&&mInT[v])throw std::runtime_error("matched atom retained in support");
    const Size owner=mBM.blossom_of(v);
    if((bool)listed[v]!=(bool)mInT[v]||((owner!=cNullItm)!=(bool)mInT[v]))
      throw std::runtime_error("support forest/BM universe mismatch");
    if(mInT[v]) {
      const Size parent=mAuxNCA.parent(v);
      if(parent>mNumVtxs||!mInT[parent]||mAuxNCA.depth(v)!=mAuxNCA.depth(parent)+1||
         mBM.mDepth[v]!=mAuxNCA.depth(v)||(parent!=mAuxRoot&&mSupportTreeOf[parent]!=mSupportTreeOf[v]))
        throw std::runtime_error("support parent/depth/component invariant");
      if(mBlsm[b].status==BS_ODD&&mBM.blossom_size(owner)!=1)
        throw std::runtime_error("inner path atom is not a singleton BM set");
      if(mBlsm[b].status==BS_EVEN&&mBM.blossom_of(mBlsm[b].base)!=owner)
        throw std::runtime_error("outer logical and support partitions differ");
    }
  }
  for(Size h=0;h<mBM.mEdges.size();++h) {
    const auto&r=mBM.mEdges[h];if(!r.inUse)continue;
    if(r.v>mNumVtxs||r.w>mNumVtxs||!mInT[r.v]||!mInT[r.w])
      throw std::runtime_error("support record references absent atom");
    Size p=r.v,steps=0;
    while(p!=r.w&&p!=cNullItm&&steps++<=mNumVtxs)p=mAuxNCA.parent(p);
    if(p!=r.w)throw std::runtime_error("support record is not an ancestor back-edge");
    const Size u=r.payload/mAdjStride,i=r.payload%mAdjStride;
    if(u>=mNumVtxs||i>=mAdj[u].size())throw std::runtime_error("support record invalid original payload");
    const bool tracked=(mBmHandleU[u][i]==h&&mBmGenU[u][i]==r.generation)||
      (mBmHandleV[u][i]==h&&mBmGenV[u][i]==r.generation);
    if(!tracked||!mD3Queued[u][i])throw std::runtime_error("support record became untracked/unqueued");
    MS_ADD(edgeChecks,1);
  }
  if(mBM.validate(false)!=0)throw std::runtime_error("support BM structural validation");
  MS_ADD(checks,1);
#else
  (void)where;
#endif
}

#ifndef MATCHBOX_GABOW_FRONTIER_AUDIT_H
#define MATCHBOX_GABOW_FRONTIER_AUDIT_H
#include <cstdio>
#include <cstring>
#include <stdexcept>
namespace Matchbox { namespace FrontierReuseAudit {
typedef unsigned long long Count;
struct State {
  Count calls,parentHits,parentFallbacks,parentLeaves,children,restored,restoredLeaves,singletons;
  Count orderedChecks,cqChecks,stableChecks,maxParentSize,maxParentCapacity,maxLiveSize,maxLiveCapacity,maxDeadCapacity,maxScratchCapacity;
  Count maxLiveMapCapacity,maxDeadMapCapacity;
  bool active;int backend;
  State() {std::memset(this,0,sizeof(*this));}
};
inline State& state() {static State s;return s;}
inline void maximum(Count& target,Count value) {if(value>target)target=value;}
struct Session {
  explicit Session(int backend) {state()=State();state().active=true;state().backend=backend;}
  ~Session() {
    State& s=state();s.active=false;
    std::printf("FRONTIER_REUSE backend=%s calls=%llu parent_hits=%llu parent_fallbacks=%llu parent_leaves=%llu children=%llu restored=%llu restored_leaves=%llu singletons=%llu order_checks=%llu cq_checks=%llu stable_checks=%llu max_parent_size=%llu max_parent_capacity=%llu max_live_size=%llu max_live_capacity=%llu max_dead_capacity=%llu max_scratch_capacity=%llu max_live_map_capacity=%llu max_dead_map_capacity=%llu\n",
      s.backend==1?"gabow":"edmonds",s.calls,s.parentHits,s.parentFallbacks,s.parentLeaves,s.children,s.restored,s.restoredLeaves,s.singletons,
      s.orderedChecks,s.cqChecks,s.stableChecks,s.maxParentSize,s.maxParentCapacity,s.maxLiveSize,s.maxLiveCapacity,s.maxDeadCapacity,s.maxScratchCapacity,s.maxLiveMapCapacity,s.maxDeadMapCapacity);
  }
};
} }
#endif

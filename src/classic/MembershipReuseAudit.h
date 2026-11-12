#ifndef MBX_MEMBERSHIP_REUSE_AUDIT_20261101_H
#define MBX_MEMBERSHIP_REUSE_AUDIT_20261101_H
#include <cstdio>
namespace Matchbox { namespace MembershipReuseAudit {
struct State {
  unsigned long long contractions=0,vertices=0,shortcuts=0,fallbacks=0;
  unsigned long long donorVertices=0,remainderVertices=0,verified=0;
  unsigned long long reusedSurfaces=0,nestedChildren=0,emptyCaches=0;
};
inline State& state(){static State s;return s;}
struct Session {
  int backend;
  explicit Session(int b):backend(b){state()=State();}
  ~Session(){const auto&s=state();std::printf(
    "MEMBERSHIP_REUSE schema=1 backend=%d contractions=%llu vertices=%llu shortcuts=%llu fallbacks=%llu donorVertices=%llu remainderVertices=%llu verified=%llu reusedSurfaces=%llu nestedChildren=%llu emptyCaches=%llu\n",
    backend,s.contractions,s.vertices,s.shortcuts,s.fallbacks,s.donorVertices,
    s.remainderVertices,s.verified,s.reusedSurfaces,s.nestedChildren,s.emptyCaches);}
};
} }
#ifdef MBX_MEMBERSHIP_REUSE_AUDIT
#define MR_ADD(field,value) do{Matchbox::MembershipReuseAudit::state().field+=(unsigned long long)(value);}while(false)
#else
#define MR_ADD(field,value) do{}while(false)
#endif
#endif

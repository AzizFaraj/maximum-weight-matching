#ifndef MBX_MONOTONE_SUPPORT_AUDIT_20261101_H
#define MBX_MONOTONE_SUPPORT_AUDIT_20261101_H
#include <cstdio>
namespace Matchbox { namespace MonotoneSupportAudit {
struct State {
  unsigned long long innerAdds=0,outerAdds=0,shrinkAdds=0,merges=0;
  unsigned long long retiredAtoms=0,retiredOwners=0,retiredTrees=0,splits=0;
  unsigned long long checks=0,edgeChecks=0,paths=0,pathAtoms=0,maxPath=0,routeLinks=0,childInspections=0;
};
inline State& state(){static State v;return v;}
struct Session {
  int backend;explicit Session(int b):backend(b){state()=State();}
  ~Session(){const auto&s=state();std::printf("MONOTONE_SUPPORT schema=1 backend=%d innerAdds=%llu outerAdds=%llu shrinkAdds=%llu merges=%llu retiredAtoms=%llu retiredOwners=%llu retiredTrees=%llu splits=%llu checks=%llu edgeChecks=%llu paths=%llu pathAtoms=%llu maxPath=%llu routeLinks=%llu childInspections=%llu\n",
    backend,s.innerAdds,s.outerAdds,s.shrinkAdds,s.merges,s.retiredAtoms,s.retiredOwners,s.retiredTrees,s.splits,s.checks,s.edgeChecks,s.paths,s.pathAtoms,s.maxPath,s.routeLinks,s.childInspections);}
};
} }
#ifdef MBX_MONOTONE_SUPPORT_AUDIT
#define MS_ADD(field,value) do{Matchbox::MonotoneSupportAudit::state().field+=(unsigned long long)(value);}while(false)
#define MS_PEAK(field,value) do{const auto msValue=(unsigned long long)(value);if(msValue>Matchbox::MonotoneSupportAudit::state().field)Matchbox::MonotoneSupportAudit::state().field=msValue;}while(false)
#else
#define MS_ADD(field,value) do{}while(false)
#define MS_PEAK(field,value) do{}while(false)
#endif
#endif

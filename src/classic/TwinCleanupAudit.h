#ifndef MATCHBOX_TWIN_CLEANUP_AUDIT_20261101_H
#define MATCHBOX_TWIN_CLEANUP_AUDIT_20261101_H
#ifdef MBX_TWIN_CLEANUP_AUDIT
#include <cstdio>
namespace Matchbox { namespace TwinCleanupAudit {
struct State {unsigned long long saved=0,batches=0,attempted=0,removed=0,stale=0,peakStaged=0,verified=0;};
inline State& state(){static State v;return v;}
struct Session {
  int backend;
  explicit Session(int b):backend(b){state()=State();}
  ~Session(){const State& s=state();
    std::printf("TWIN_CLEANUP schema=1 backend=%d saved=%llu batches=%llu attempted=%llu removed=%llu stale=%llu peakStaged=%llu verified=%llu\n",
      backend,s.saved,s.batches,s.attempted,s.removed,s.stale,s.peakStaged,s.verified);}
};
} }
#define TC_ADD(field,value) do{Matchbox::TwinCleanupAudit::state().field+=(unsigned long long)(value);}while(false)
#define TC_PEAK(value) do{const auto tcValue=(unsigned long long)(value);if(tcValue>Matchbox::TwinCleanupAudit::state().peakStaged)Matchbox::TwinCleanupAudit::state().peakStaged=tcValue;}while(false)
#else
#define TC_ADD(field,value) do{}while(false)
#define TC_PEAK(value) do{}while(false)
#endif
#endif

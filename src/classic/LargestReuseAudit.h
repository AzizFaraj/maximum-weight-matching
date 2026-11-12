#ifndef MBX_LARGEST_REUSE_AUDIT_20261030_H
#define MBX_LARGEST_REUSE_AUDIT_20261030_H
#include <cstdio>
#include <stdexcept>
#include <vector>
#include <limits>
namespace Matchbox { namespace LargestReuseAudit {
inline void require(bool ok, const char* why) {
  if (!ok) throw std::runtime_error(why);
}
#define LR_FIELDS(X) X(calls) X(parentLeaves) X(largestLeaves) X(smallerCopied) \
  X(compactedCopied) X(largestMoved) X(bmPartitions) X(bmParentLeaves) \
  X(bmSmallerLeaves) X(bmLargestLeaves) X(preChecks) X(postChecks)
struct State {
#define LR_DECLARE(name) unsigned long long name;
  LR_FIELDS(LR_DECLARE)
#undef LR_DECLARE
  State() {
#define LR_ZERO(name) name = 0;
    LR_FIELDS(LR_ZERO)
#undef LR_ZERO
  }
};
inline State& state() { static State value; return value; }
inline void add(unsigned long long& target, unsigned long long value) {
  require(value <= (std::numeric_limits<unsigned long long>::max)() - target, "largest audit overflow");
  target += value;
}
struct Session {
  int backend;
  explicit Session(int b) : backend(b) { state()=State(); }
  ~Session() {
    std::printf("LARGEST_REUSE schema=1 backend=%d", backend);
#define LR_EMIT(name) std::printf(" " #name "=%llu", state().name);
    LR_FIELDS(LR_EMIT)
#undef LR_EMIT
    std::printf("\n");
  }
};
} }
#ifdef MBX_LARGEST_REUSE_AUDIT
#define LR_ADD(name, amount) Matchbox::LargestReuseAudit::add(Matchbox::LargestReuseAudit::state().name,(amount))
#else
#define LR_ADD(name, amount) ((void)0)
#endif
#endif

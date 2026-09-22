#ifndef SCALING_BUCKETPQ_H
#define SCALING_BUCKETPQ_H

// Bucket priority queue for the S5 event-driven Edmonds search (paper §5.2,
// "BucketSearch"). Events are scheduled at integer times (= the dual-adjustment
// count at which they fire). bucket[i] = events scheduled for time i. All ops are
// O(1) amortized; advancing the clock scans empty buckets in O(t_max) total over
// a search. When the max time is unbounded, a real PQSearch would swap this for a
// general integer PQ -- but for BucketSearch (bounded dual adjustments per scale)
// the bucket array is exact and O(1).
//
// This is a standalone, independently-tested component so the search engine is
// built on a verified PQ. Events carry an opaque payload (the edge / blossom id);
// the search interprets it. Stale events (whose payload is no longer relevant) are
// filtered by the search when popped (lazy deletion), as the paper does.

#include <vector>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <utility>

namespace Matchbox {

class BucketPQ {
public:
  typedef long long LL;
  struct Event { LL time; long long a; long long b; int kind; long long owner; }; // owner: source outer node
  // Bounded searches dominate Hybrid's queue traffic. Keep their retained
  // records compact while expanding to the wide public Event at pop time.
  // UINT32_MAX represents the existing -1 sentinel; all other identifiers
  // remain exact through 2^32-2 (far beyond an in-memory million-node graph).
  struct StoredEvent {
    LL time;
    std::uint32_t a, b, owner;
    std::uint8_t kind;
    StoredEvent(const Event& event) : time(event.time),
      a(pack(event.a)), b(pack(event.b)), owner(pack(event.owner)),
      kind((std::uint8_t)event.kind) {}
    static std::uint32_t pack(long long value) {
      return value < 0 ? std::numeric_limits<std::uint32_t>::max() :
                         (std::uint32_t)value;
    }
    static long long unpack(std::uint32_t value) {
      return value == std::numeric_limits<std::uint32_t>::max() ? -1LL :
                                                                  (long long)value;
    }
    void expand(Event* event) const {
      event->time = time; event->a = unpack(a); event->b = unpack(b);
      event->kind = (int)kind; event->owner = unpack(owner);
    }
  };
  static_assert(sizeof(StoredEvent) <= 24,
                "bounded event payload unexpectedly lost its compact layout");

  BucketPQ() : mMinOcc(MAXT), mMaxOcc(-1), mNow(0), mCount(0), mGen(0),
               mRetainedEventSlots(0) {}

  // capacity = max event time + 1 (exclusive upper bound on scheduled times).
  void init(LL capacity) {
    if (capacity < 1) capacity = 1;
    mBuckets.assign((size_t)capacity, std::vector<StoredEvent>());
    mHeads.assign((size_t)capacity, 0);
    mTouched.clear(); mMinOcc = MAXT; mMaxOcc = -1;
    mNow = 0; mCount = 0; mRetainedEventSlots = 0;
  }

  LL now() const { return mNow; }
  // Re-sync the scan cursor DOWN to the caller's logical clock. advanceToNext PEEKS ahead to the
  // next occupied bucket (mNow may exceed the true search time); an insert made after that peek
  // but at the still-current logical time would otherwise be clamped FORWARD by schedule() and
  // fire late (the persistent-boundary re-offer bug: events re-inserted at t=93 popped at t=106
  // where strict eligibility discarded them). Rewinding the cursor is always safe: it only makes
  // advanceToNext re-scan a few buckets.
  void syncNow(LL t) { if (t >= 0 && t < mNow) mNow = t; }
  bool empty() const { return mCount == 0; }
  size_t size() const { return mCount; }

  // schedule an event to fire at absolute time t (>= mNow). t must be < capacity.
  bool schedule(LL t, long long a, long long b, int kind, long long owner = -1) {
    const unsigned long long maxId =
      (unsigned long long)std::numeric_limits<std::uint32_t>::max() - 1ULL;
    if (a < -1 || b < -1 || owner < -1 || kind < 0 || kind > 255 ||
        (a >= 0 && (unsigned long long)a > maxId) ||
        (b >= 0 && (unsigned long long)b > maxId) ||
        (owner >= 0 && (unsigned long long)owner > maxId)) return false;
    if (t < mNow) t = mNow;                 // never schedule in the past
    Event e; e.time = t; e.a = a; e.b = b; e.kind = kind; e.owner = owner;
    StoredEvent stored(e);
    // OVERFLOW GUARD: the bucket array is O(max scheduled time). With carried Property-2 duals
    // (y<-2y'+3 each scale) slacks can blow up ~2^scale, so a far-future event would try to
    // allocate O(2^scale) buckets (measured: 32GB OOM in the bridge FVR). Cap the bucket window;
    // events beyond mWindowCap go into an ORDERED overflow set (time-keyed), pulled into buckets
    // lazily as the clock advances. Keeps the common (near) case O(1) and never over-allocates.
    // This is the "general integer PQ for unbounded t" the header notes -- localized to the tail.
    if ((size_t)t >= mWindowCap) {
      mOverflow.insert(std::make_pair(t, stored));
      if (t < mMinOcc) mMinOcc = t;
      if (t > mMaxOcc) mMaxOcc = t;
      ++mCount;
      return true;
    }
    if ((size_t)t >= mBuckets.size()) growTo(t);
    // FIFO within a dual level.  Figure-6 permits arbitrary tied-event order, but preserving
    // insertion order also matches the scan engine's deterministic adjacency traversal.  mHeads
    // makes pop O(1); when a bucket drains, popAtNow clears it while retaining its capacity.
    if (mHeads[(size_t)t] >= mBuckets[(size_t)t].size()) {
      mBuckets[(size_t)t].clear(); mHeads[(size_t)t] = 0;
      mTouched.push_back(t);                       // O(1) record (dups harmless)
    }
    size_t oldCapacity = mBuckets[(size_t)t].capacity();
    mBuckets[(size_t)t].push_back(stored);
    mRetainedEventSlots += mBuckets[(size_t)t].capacity() - oldCapacity;
    if (t < mMinOcc) mMinOcc = t;           // min-occupied hint -> O(1) skip of empty prefix
    if (t > mMaxOcc) mMaxOcc = t;
    ++mCount;
    return true;
  }

  // advance the clock to the next non-empty bucket and return whether one exists.
  // Start from max(mNow, mMinOcc) so we skip the leading run of empty buckets in O(1)
  // (the common case: mNow lags far behind the first real event after a relabel/jump).
  // The residual scan only covers gaps BETWEEN real events, which are sparse.
  bool advanceToNext() {
#ifdef S5DBG
    extern long long gBPQadvanceCalls; ++gBPQadvanceCalls;
#endif
    if (mCount == 0) return false;
    if (mNow < mMinOcc) mNow = mMinOcc;     // jump over the empty prefix
    // scan bucket window for the next occupied bucket (bounded by mWindowCap)
    while ((size_t)mNow <= (size_t)mMaxOcc && (size_t)mNow < mBuckets.size()
           && (size_t)mNow < mWindowCap
           && mHeads[(size_t)mNow] >= mBuckets[(size_t)mNow].size()) {
#ifdef S5DBG
      { extern long long gBPQemptyScan; ++gBPQemptyScan; }
#endif
      ++mNow;
    }
    bool inWindow = (size_t)mNow < mBuckets.size() && (size_t)mNow < mWindowCap
                    && mHeads[(size_t)mNow] < mBuckets[(size_t)mNow].size();
    if (inWindow) return true;
    // window exhausted -> the next event (if any) is in the overflow set.
    if (!mOverflow.empty()) { mNow = mOverflow.begin()->first; return true; }
    return false;
  }

  // pop one event at the current time (call advanceToNext first). Returns false if
  // the current time has no event.
  bool popAtNow(Event* out) {
    if ((size_t)mNow < mBuckets.size() && (size_t)mNow < mWindowCap
        && mHeads[(size_t)mNow] < mBuckets[(size_t)mNow].size()) {
      std::vector<StoredEvent>& bucket = mBuckets[(size_t)mNow];
      size_t& head = mHeads[(size_t)mNow];
      bucket[head++].expand(out);
      if (head == bucket.size()) { bucket.clear(); head = 0; }
      --mCount;
      if (mCount == 0) { mMinOcc = MAXT; mMaxOcc = -1; }
      return true;
    }
    // overflow event at the current time?
    if (!mOverflow.empty() && mOverflow.begin()->first == mNow) {
      mOverflow.begin()->second.expand(out);
      mOverflow.erase(mOverflow.begin());
      --mCount;
      if (mCount == 0) { mMinOcc = MAXT; mMaxOcc = -1; }
      return true;
    }
    return false;
  }

  void clear() { ++mGen; mBuckets.clear(); mHeads.clear(); mTouched.clear(); mOverflow.clear(); mMinOcc = MAXT; mMaxOcc = -1; mNow = 0; mCount = 0; mRetainedEventSlots = 0; }

  // reset for reuse WITHOUT freeing the bucket-array storage (the per-phase clear()+
  // init() was reallocating the whole array every phase -- a hot-path alloc). Empties
  // each bucket up to the high-water mark in place (keeps each bucket's capacity too).
  // Only touches buckets that were occupied (via mOccupied) so this is O(t), not O(maxW).
  // GENERATION: bumps on every reset()/clear() -- lets a caller detect whether events it scheduled
  // are still live (same gen) or were wiped (gen changed). Used by PS5_INCR_EVENTS to skip re-emitting
  // a big blossom's still-live events on a marker-pop that had no reset since the schedule.
  long long generation() const { return mGen; }
  void reset() {
    ++mGen;
    // Reusing modest bucket capacities avoids hot-path allocation. Across many
    // independent searches, however, different time buckets can each retain a
    // large historical peak; that accumulation reached tens of GiB at n=100k.
    // Once retained Event storage exceeds this ceiling, destroy the inner
    // vectors and return to a one-bucket queue.
    if (mRetainedEventSlots > MAX_RETAINED_EVENT_SLOTS) {
      mBuckets.assign(1, std::vector<StoredEvent>());
      mHeads.assign(1, 0);
      mTouched.clear(); mOverflow.clear(); mMinOcc = MAXT; mMaxOcc = -1;
      mNow = 0; mCount = 0; mRetainedEventSlots = 0;
      return;
    }
    // Clear ONLY the buckets that were touched (O(#schedules), not O(maxW) span).
    for (size_t i = 0; i < mTouched.size(); ++i) {
      size_t t = (size_t)mTouched[i];
      if (t < mBuckets.size()) { mBuckets[t].clear(); mHeads[t] = 0; }
    }
    mTouched.clear(); mOverflow.clear(); mMinOcc = MAXT; mMaxOcc = -1; mNow = 0; mCount = 0;
  }

private:
  void growTo(LL t) {
    size_t need = (size_t)t + 1;
    if (need > mBuckets.size()) { mBuckets.resize(need); mHeads.resize(need, 0); }
  }
  static const LL MAXT = (LL)1 << 62;
  // bucket window: times in [0,mWindowCap) use the O(1) bucket array; times >= mWindowCap go into
  // mOverflow (ordered, O(log)). Caps bucket memory at ~mWindowCap*sizeof(vector) regardless of how
  // large a scheduled time gets (the carried-dual slack blowup). 1<<20 buckets ~= a few MB.
  static const size_t mWindowCap = (size_t)1 << 20;
  // About 96 MiB with the compact 24-byte retained-event layout.
  static const size_t MAX_RETAINED_EVENT_SLOTS = (size_t)4 << 20;
  std::vector<std::vector<StoredEvent> > mBuckets;
  std::vector<size_t> mHeads;       // next FIFO event in each bucket
  std::vector<LL> mTouched;   // times push_back'd to since last reset (for O(t) reset)
  struct OvCmp { bool operator()(const std::pair<LL, StoredEvent>& x, const std::pair<LL, StoredEvent>& y) const { return x.first < y.first; } };
  std::multiset<std::pair<LL, StoredEvent>, OvCmp> mOverflow;  // far-future events (time >= mWindowCap), time-ordered
  LL mMinOcc, mMaxOcc;        // occupied-time hints (skip empty prefix; bound the scan)
  LL mNow;
  size_t mCount;
  long long mGen;             // bumps on reset()/clear(); lets callers detect wiped events
  size_t mRetainedEventSlots;
};

} // namespace Matchbox

#endif // SCALING_BUCKETPQ_H

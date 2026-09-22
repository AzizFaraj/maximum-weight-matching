#ifndef MATCHBOX_BINARY_EVENT_HEAP_H
#define MATCHBOX_BINARY_EVENT_HEAP_H

#include <vector>
#include <cstddef>

namespace Matchbox {

// Unbounded event queue used by practical PQSearch. Equal-time events retain
// insertion order so scan and event backends have deterministic tie behavior.
class BinaryEventHeap {
public:
  typedef long long Time;
  struct Event {
    Time time;
    long long a, b;
    int kind;
    long long owner;
    unsigned long long serial;
    Event() : time(0), a(0), b(0), kind(0), owner(-1), serial(0) {}
  };

  BinaryEventHeap() : mSerial(0), mBatching(false) {}
  void clear() { mHeap.clear(); mSerial = 0; mBatching = false; }
  bool empty() const { return mHeap.empty(); }
  std::size_t size() const { return mHeap.size(); }

  // Initial PQSearch seeding knows the complete O(m) event set. Append it
  // unsifted and build the heap once in linear time; later dynamic events use
  // ordinary logarithmic insertion.
  void beginBatch() { mBatching = true; }
  void endBatch() {
    if (!mBatching) return;
    for (std::size_t i = mHeap.size() / 2; i > 0; --i) siftDown(i - 1);
    mBatching = false;
  }

  void schedule(Time time, long long a, long long b, int kind,
                long long owner = -1) {
    Event event; event.time = time; event.a = a; event.b = b;
    event.kind = kind; event.owner = owner; event.serial = mSerial++;
    mHeap.push_back(event);
    if (mBatching) return;
    std::size_t i = mHeap.size() - 1;
    while (i > 0) {
      std::size_t p = (i - 1) / 2;
      if (!less(mHeap[i], mHeap[p])) break;
      Event tmp = mHeap[i]; mHeap[i] = mHeap[p]; mHeap[p] = tmp;
      i = p;
    }
  }

  bool pop(Event* out) {
    if (!out || mHeap.empty()) return false;
    *out = mHeap[0];
    Event tail = mHeap.back(); mHeap.pop_back();
    if (mHeap.empty()) return true;
    mHeap[0] = tail;
    siftDown(0);
    return true;
  }

private:
  void siftDown(std::size_t i) {
    for (;;) {
      std::size_t left = 2 * i + 1, right = left + 1, best = i;
      if (left < mHeap.size() && less(mHeap[left], mHeap[best])) best = left;
      if (right < mHeap.size() && less(mHeap[right], mHeap[best])) best = right;
      if (best == i) break;
      Event tmp = mHeap[i]; mHeap[i] = mHeap[best]; mHeap[best] = tmp;
      i = best;
    }
  }
  static bool less(const Event& x, const Event& y) {
    return x.time < y.time || (x.time == y.time && x.serial < y.serial);
  }
  std::vector<Event> mHeap;
  unsigned long long mSerial;
  bool mBatching;
};

} // namespace Matchbox

#endif

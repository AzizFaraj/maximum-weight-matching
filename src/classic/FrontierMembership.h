#ifndef MBX_FRONTIER_MEMBERSHIP_20261029_H
#define MBX_FRONTIER_MEMBERSHIP_20261029_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <class Size> class FrontierMembership;

// Non-owning read-only membership. Views never extend allocation lifetimes;
// discard them before their owner mutates, moves, releases or is destroyed.
template <class Size>
class FrontierMembershipView {
public:
  typedef Size value_type;
  typedef std::size_t size_type;
  typedef const Size* const_iterator;

  FrontierMembershipView() noexcept : data_(nullptr), size_(0) {}
  FrontierMembershipView(const Size* data, size_type size)
      : data_(data), size_(size) {
    if (size != 0 && data == nullptr)
      throw std::invalid_argument("null nonempty membership view");
  }
  FrontierMembershipView(const std::vector<Size>& values) noexcept
      : data_(values.data()), size_(values.size()) {}
  FrontierMembershipView(const std::vector<Size>* values) noexcept
      : data_(values ? values->data() : nullptr), size_(values ? values->size() : 0) {}
  FrontierMembershipView(const FrontierMembership<Size>& values) noexcept;
  FrontierMembershipView(const FrontierMembership<Size>* values) noexcept;
  const FrontierMembershipView* operator->() const noexcept { return this; }
  const FrontierMembershipView& operator*() const noexcept { return *this; }

  size_type size() const noexcept { return size_; }
  bool empty() const noexcept { return size_ == 0; }
  const Size* data() const noexcept { return data_; }
  const_iterator begin() const noexcept { return data_; }
  const_iterator end() const noexcept {
    return size_ == 0 ? data_ : data_ + size_;
  }
  const Size& operator[](size_type index) const noexcept {
    assert(index < size_);
    return data_[index];
  }
  const Size& at(size_type index) const {
    if (index >= size_) throw std::out_of_range("membership view index");
    return data_[index];
  }
  const Size& back() const noexcept {
    assert(!empty());
    return data_[size_ - 1];
  }

private:
  const Size* data_;
  size_type size_;
};

// One uniquely owned vector plus an active-start offset, with NO cached length.
// The active range is ALWAYS the suffix [offset_, backing_.size()). takeRange
// destroys any unselected suffix before publication; for the engine's integral
// Size this shrinking resize does not allocate, move or copy surviving values.
// Smaller children must be copied before a largest-child takeRange consumes the
// parent. No shared ownership can pin a historical ancestor allocation.
//
// capacity() reports the ENTIRE allocation, including inactive prefix. A range
// transfer, deep copy and clear enforce retentionLimit(size()). Explicit reserve
// can exceed it during construction, as before. Direct push/append deliberately
// delegate growth to vector; callers growing an arbitrary prefixed slice can
// call compactIfNeeded() at their publication boundary. Engine contraction uses
// fresh offset-zero parents, keeping its hot push identical to vector push.
template <class Size>
class FrontierMembership {
public:
  typedef Size value_type;
  typedef std::size_t size_type;
  typedef Size* iterator;
  typedef const Size* const_iterator;
  typedef FrontierMembershipView<Size> View;

  FrontierMembership() noexcept : offset_(0) {}
  explicit FrontierMembership(const std::vector<Size>& values)
      : backing_(values), offset_(0) {}
  FrontierMembership(std::initializer_list<Size> values)
      : backing_(values), offset_(0) {}
  template <class InputIt,
            typename std::enable_if<!std::is_integral<InputIt>::value,
                                    int>::type = 0>
  FrontierMembership(InputIt first, InputIt last) : offset_(0) {
    assign(first, last);
  }

  // Copies contain only the logical active suffix, never its inactive prefix.
  FrontierMembership(const FrontierMembership& other) : offset_(0) {
    if (!other.empty()) backing_.assign(other.begin(), other.end());
  }
  FrontierMembership& operator=(const FrontierMembership& other) {
    if (this != &other) {
      FrontierMembership copy(other);
      swap(copy);
    }
    return *this;
  }
  FrontierMembership(FrontierMembership&& other) noexcept
      : backing_(std::move(other.backing_)), offset_(other.offset_) {
    other.release();
  }
  FrontierMembership& operator=(FrontierMembership&& other) noexcept {
    if (this != &other) {
      FrontierMembership moved(std::move(other));
      swap(moved);
    }
    return *this;
  }

  size_type size() const noexcept { return backing_.size() - offset_; }
  bool empty() const noexcept { return backing_.size() == offset_; }
  size_type capacity() const noexcept { return backing_.capacity(); }
  size_type backingOffset() const noexcept { return offset_; }
  static size_type retentionLimit(size_type liveSize) noexcept {
    const size_type maximum = (std::numeric_limits<size_type>::max)();
    return liveSize > (maximum - size_type(4096)) / size_type(4)
               ? maximum : size_type(4) * liveSize + size_type(4096);
  }
  bool storageWithinBound() const noexcept {
    return capacity() <= retentionLimit(size());
  }

  Size* data() noexcept {
    return offset_ == 0 ? backing_.data() : backing_.data() + offset_;
  }
  const Size* data() const noexcept {
    return offset_ == 0 ? backing_.data() : backing_.data() + offset_;
  }
  iterator begin() noexcept { return data(); }
  const_iterator begin() const noexcept { return data(); }
  const_iterator cbegin() const noexcept { return data(); }
  iterator end() noexcept {
    return backing_.empty() ? backing_.data() : backing_.data() + backing_.size();
  }
  const_iterator end() const noexcept {
    return backing_.empty() ? backing_.data() : backing_.data() + backing_.size();
  }
  const_iterator cend() const noexcept { return end(); }
  Size& operator[](size_type index) noexcept {
    assert(index < size());
    return backing_[offset_ + index];
  }
  const Size& operator[](size_type index) const noexcept {
    assert(index < size());
    return backing_[offset_ + index];
  }
  Size& at(size_type index) {
    if (index >= size()) throw std::out_of_range("membership index");
    return (*this)[index];
  }
  const Size& at(size_type index) const {
    if (index >= size()) throw std::out_of_range("membership index");
    return (*this)[index];
  }
  Size& back() noexcept {
    assert(!empty());
    return backing_.back();
  }
  const Size& back() const noexcept {
    assert(!empty());
    return backing_.back();
  }

  void clear() noexcept {
    backing_.clear();
    offset_ = 0;
    if (!storageWithinBound()) release();
  }
  void release() noexcept {
    std::vector<Size>().swap(backing_);
    offset_ = 0;
  }
  void swap(FrontierMembership& other) noexcept {
    backing_.swap(other.backing_);
    std::swap(offset_, other.offset_);
  }
  // requested is active capacity, while capacity() reports full allocation.
  void reserve(size_type requested) {
    if (requested <= backing_.capacity() - offset_) return;
    normalizePrefix();
    backing_.reserve(requested);
  }
  // No logical-size bookkeeping, offset test, inactive-tail branch or reserve.
  // By-value input also protects push_back(membership[i]) during reallocation.
  void push_back(Size value) { backing_.push_back(std::move(value)); }

  template <class InputIt,
            typename std::enable_if<!std::is_integral<InputIt>::value,
                                    int>::type = 0>
  void assign(InputIt first, InputIt last) {
    if (first == last) {
      clear();
      return;
    }
    std::vector<Size> replacement(first, last);
    backing_.swap(replacement);
    offset_ = 0;
  }
  void assign(size_type count, const Size& value) {
    std::vector<Size> replacement(count, value);
    backing_.swap(replacement);
    offset_ = 0;
  }
  void assign(const std::vector<Size>& values) {
    assign(values.begin(), values.end());
  }
  void assign(View values) { assign(values.begin(), values.end()); }

  void append(const std::vector<Size>& values) { append(View(values)); }
  void append(const FrontierMembership& values) { append(View(values)); }
  void append(View values) {
    if (values.empty()) return;
    // Physical vector size includes the prefix, so check physical growth.
    if (values.size() > backing_.max_size() - backing_.size())
      throw std::length_error("membership append overflow");
    if (aliasesBacking(values.data())) {
      std::vector<Size> copy(values.begin(), values.end());
      append(View(copy));
      return;
    }
    // An active suffix ends at the physical vector end even with an offset.
    backing_.insert(backing_.end(), values.begin(), values.end());
  }

  // Invalid/self ranges leave both objects and copiedElements unchanged.
  // Success consumes the whole parent; a zero range releases both objects.
  // copiedElements counts only range compaction, not the separately copied
  // smaller children. Engine Size has trivial, nonthrowing destruction, making
  // the zero-copy branch's shrinking resize nonallocating and nonthrowing.
  bool takeRange(FrontierMembership& parent, size_type offset, size_type count,
                 size_type* copiedElements = nullptr) {
    if (this == &parent || offset > parent.size() ||
        count > parent.size() - offset) return false;
    size_type copied = 0;
    if (count == 0) {
      release();
      parent.release();
    } else if (parent.capacity() > retentionLimit(count)) {
      // Allocate before mutation: failure preserves target and parent.
      FrontierMembership compact(parent.begin() + offset,
                                  parent.begin() + offset + count);
      swap(compact);
      parent.release();
      copied = count;
    } else {
      // Bounds above prove both sums fit and selectedEnd <= old vector size.
      const size_type selectedOffset = parent.offset_ + offset;
      const size_type selectedEnd = selectedOffset + count;
      *this = std::move(parent);
      backing_.resize(selectedEnd);
      offset_ = selectedOffset;
    }
    if (copiedElements != nullptr) *copiedElements = copied;
    return true;
  }

  size_type compactIfNeeded() {
    if (storageWithinBound()) return 0;
    const size_type copied = size();
    FrontierMembership compact(begin(), end());
    swap(compact);
    return copied;
  }

private:
  void normalizePrefix() {
    if (offset_ == 0) return;
    const size_type active = size();
    std::move(backing_.begin() + offset_, backing_.end(), backing_.begin());
    backing_.resize(active);
    offset_ = 0;
  }
  bool aliasesBacking(const Size* pointer) const noexcept {
    if (backing_.empty()) return false;
    const Size* first = backing_.data();
    const Size* last = first + backing_.size();
    const std::less<const Size*> less;
    return !less(pointer, first) && less(pointer, last);
  }

  std::vector<Size> backing_;
  size_type offset_;
};

template <class Size>
FrontierMembershipView<Size>::FrontierMembershipView(
    const FrontierMembership<Size>& values) noexcept
    : data_(values.data()), size_(values.size()) {}

template <class Size>
FrontierMembershipView<Size>::FrontierMembershipView(
    const FrontierMembership<Size>* values) noexcept
    : data_(values ? values->data() : nullptr), size_(values ? values->size() : 0) {}

template <class Size>
void swap(FrontierMembership<Size>& left,
          FrontierMembership<Size>& right) noexcept {
  left.swap(right);
}

#endif

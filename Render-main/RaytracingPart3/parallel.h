#pragma once
#include "geo.h"
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include "efloat.h"
#include "check.h"


// Parallel Declarations
class AtomicFloat {
public:
	// AtomicFloat Public Methods
	explicit AtomicFloat(Float v = 0) { bits = FloatToBits(v); }
	operator Float() const { return BitsToFloat(bits); }
	Float operator=(Float v) {
		bits = FloatToBits(v);
		return v;
	}
	void Add(Float v) {
#ifdef PBRT_FLOAT_AS_DOUBLE
		uint64_t oldBits = bits, newBits;
#else
		uint32_t oldBits = bits, newBits;
#endif
		do {
			newBits = FloatToBits(BitsToFloat(oldBits) + v);
		} while (!bits.compare_exchange_weak(oldBits, newBits));
	}

private:
	// AtomicFloat Private Data
#ifdef PBRT_FLOAT_AS_DOUBLE
	std::atomic<uint64_t> bits;
#else
	std::atomic<uint32_t> bits;
#endif
};

// Simple one-use barrier; ensures that multiple threads all reach a
// particular point of execution before allowing any of them to proceed
// past it.
//
// Note: this should be heap allocated and managed with a shared_ptr, where
// all threads that use it are passed the shared_ptr. This ensures that
// memory for the Barrier won't be freed until all threads have
// successfully cleared it.
class Barrier {
public:
	Barrier(int count) : count(count) { check::CHECK_GT(count, 0); }
	~Barrier() { check::CHECK_EQ(count, 0); }
	void Wait();

private:
	std::mutex mutex;
	std::condition_variable cv;
	int count;
};

void ParallelFor(std::function<void(int64_t)> func, int64_t count,
	int chunkSize = 1);
extern thread_local int ThreadIndex;
void ParallelFor2D(std::function<void(Point2i)> func, const Point2i& count);
int MaxThreadIndex();
int NumSystemCores();

void ParallelInit();
void ParallelCleanup();
void MergeWorkerThreadStats();
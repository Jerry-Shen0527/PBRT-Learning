#include "stats.h"
thread_local uint64_t ProfilerState;
std::vector<std::function<void(StatsAccumulator&)>>* StatRegisterer::funcs;
static StatsAccumulator statsAccumulator;
void ProfilerWorkerThreadInit() {
#ifdef PBRT_HAVE_ITIMER
    // The per-thread initialization in the worker threads has to happen
    // *before* the profiling signal handler is installed.
    CHECK(!profilerRunning || profilerSuspendCount > 0);

    // ProfilerState is a thread-local variable that is accessed in the
    // profiler signal handler. It's important to access it here, which
    // causes the dynamic memory allocation for the thread-local storage to
    // happen now, rather than in the signal handler, where this isn't
    // allowed.
    ProfilerState = ProfToBits(Prof::SceneConstruction);
#endif  // PBRT_HAVE_ITIMER
}
void ReportThreadStats() {
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    StatRegisterer::CallCallbacks(statsAccumulator);
}
void StatRegisterer::CallCallbacks(StatsAccumulator& accum) {
    for (auto func : *funcs) func(accum);
}
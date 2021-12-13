#pragma once
#include "efloat.h"
#include <map>
#include <cfloat>
#include <chrono>
#include <string>
#include <functional>
#include <mutex>
#include <vector>
//#ifndef ProfilerState
//thread_local uint64_t ProfilerState;
//#endif // !1
using namespace std;
extern thread_local uint64_t ProfilerState;

enum class Prof {
    SceneConstruction,
    AccelConstruction,
    TextureLoading,
    MIPMapCreation,

    IntegratorRender,
    SamplerIntegratorLi,
    SPPMCameraPass,
    SPPMGridConstruction,
    SPPMPhotonPass,
    SPPMStatsUpdate,
    BDPTGenerateSubpath,
    BDPTConnectSubpaths,
    LightDistribLookup,
    LightDistribSpinWait,
    LightDistribCreation,
    DirectLighting,
    BSDFEvaluation,
    BSDFSampling,
    BSDFPdf,
    BSSRDFEvaluation,
    BSSRDFSampling,
    PhaseFuncEvaluation,
    PhaseFuncSampling,
    AccelIntersect,
    AccelIntersectP,
    LightSample,
    LightPdf,
    MediumSample,
    MediumTr,
    TriIntersect,
    TriIntersectP,
    CurveIntersect,
    CurveIntersectP,
    ShapeIntersect,
    ShapeIntersectP,
    ComputeScatteringFuncs,
    GenerateCameraRay,
    MergeFilmTile,
    SplatFilm,
    AddFilmSample,
    StartPixel,
    GetSample,
    TexFiltTrilerp,
    TexFiltEWA,
    TexFiltPtex,
    NumProfCategories
};
inline uint64_t ProfToBits(Prof p) { return 1ull << (int)p; }
class ProfilePhase {
public:
    // ProfilePhase Public Methods
    ProfilePhase(Prof p) {
        categoryBit = ProfToBits(p);
        reset = (ProfilerState & categoryBit) == 0;
        ProfilerState |= categoryBit;
    }
    ~ProfilePhase() {
        if (reset) ProfilerState &= ~categoryBit;
    }
    ProfilePhase(const ProfilePhase&) = delete;
    ProfilePhase& operator=(const ProfilePhase&) = delete;

private:
    // ProfilePhase Private Data
    bool reset;
    uint64_t categoryBit;
};
class StatsAccumulator;
class StatRegisterer {
public:
    // StatRegisterer Public Methods
    StatRegisterer(std::function<void(StatsAccumulator&)> func) {
        static std::mutex mutex;
        std::lock_guard<std::mutex> lock(mutex);
        if (!funcs)
            funcs = new std::vector<std::function<void(StatsAccumulator&)>>;
        funcs->push_back(func);
    }
    static void CallCallbacks(StatsAccumulator& accum);

private:
    // StatRegisterer Private Data
    static std::vector<std::function<void(StatsAccumulator&)>>* funcs;
};

class StatsAccumulator {
public:
    // StatsAccumulator Public Methods
    void ReportCounter(const std::string& name, int64_t val) {
        counters[name] += val;
    }
    void ReportMemoryCounter(const std::string& name, int64_t val) {
        memoryCounters[name] += val;
    }
    void ReportIntDistribution(const std::string& name, int64_t sum,
        int64_t count, int64_t min, int64_t max) {
        intDistributionSums[name] += sum;
        intDistributionCounts[name] += count;
        if (intDistributionMins.find(name) == intDistributionMins.end())
            intDistributionMins[name] = min;
        else
            intDistributionMins[name] =
            std::min(intDistributionMins[name], min);
        if (intDistributionMaxs.find(name) == intDistributionMaxs.end())
            intDistributionMaxs[name] = max;
        else
            intDistributionMaxs[name] =
            std::max(intDistributionMaxs[name], max);
    }
    void ReportFloatDistribution(const std::string& name, double sum,
        int64_t count, double min, double max) {
        floatDistributionSums[name] += sum;
        floatDistributionCounts[name] += count;
        if (floatDistributionMins.find(name) == floatDistributionMins.end())
            floatDistributionMins[name] = min;
        else
            floatDistributionMins[name] =
            std::min(floatDistributionMins[name], min);
        if (floatDistributionMaxs.find(name) == floatDistributionMaxs.end())
            floatDistributionMaxs[name] = max;
        else
            floatDistributionMaxs[name] =
            std::max(floatDistributionMaxs[name], max);
    }
    void ReportPercentage(const std::string& name, int64_t num, int64_t denom) {
        percentages[name].first += num;
        percentages[name].second += denom;
    }
    void ReportRatio(const std::string& name, int64_t num, int64_t denom) {
        ratios[name].first += num;
        ratios[name].second += denom;
    }

    void Print(FILE* file);
    void Clear();

private:
    // StatsAccumulator Private Data
    std::map<std::string, int64_t> counters;
    std::map<std::string, int64_t> memoryCounters;
    std::map<std::string, int64_t> intDistributionSums;
    std::map<std::string, int64_t> intDistributionCounts;
    std::map<std::string, int64_t> intDistributionMins;
    std::map<std::string, int64_t> intDistributionMaxs;
    std::map<std::string, double> floatDistributionSums;
    std::map<std::string, int64_t> floatDistributionCounts;
    std::map<std::string, double> floatDistributionMins;
    std::map<std::string, double> floatDistributionMaxs;
    std::map<std::string, std::pair<int64_t, int64_t>> percentages;
    std::map<std::string, std::pair<int64_t, int64_t>> ratios;
};

//static StatsAccumulator statsAccumulator;
void ProfilerWorkerThreadInit();
void ReportThreadStats(); 
//#ifndef ProfilerState
//extern thread_local uint64_t ProfilerState;
//#endif // ProfilerState


inline uint64_t CurrentProfilerState() { return ProfilerState; }

#define STAT_COUNTER(title, var)                           \
    static thread_local int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportCounter(title, var);                   \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)
#define STAT_MEMORY_COUNTER(title, var)                    \
    static thread_local int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportMemoryCounter(title, var);             \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)
#define STAT_RATIO(title, numVar, denomVar)                   \
    static thread_local int64_t numVar, denomVar;        \
    static void STATS_FUNC##numVar(StatsAccumulator &accum) { \
        accum.ReportRatio(title, numVar, denomVar);           \
        numVar = denomVar = 0;                                \
    }                                                         \
    static StatRegisterer STATS_REG##numVar(STATS_FUNC##numVar)
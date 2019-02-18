#pragma once

#include <map>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstddef>             // for size_t
#include <cmath>               // for fabs

#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>

std::vector<double> sg_smooth(const std::vector<double> &v, const int w, const int deg);
std::vector<double> sg_derivative(const std::vector<double> &v, const int w,
                                  const int deg, const double h = 1.0);

class SGSmoothing {
public:
    static void Smooth(size_t nsamples, double* in, double* out, const int w, const int deg);
    static void Derivative(size_t nsamples, double* in, double* out, const int w, const int deg, const double h = 1.0);
};

enum CFDParams {
    kADCNBits,
    kADCDynamicRange,
    kADCOffset,
    kTimeSamplingInterval,
    kNSamplingPoints,
    kCFDThreshold,
    kRawHitFinderThresholdInNoiseSigma,
    kShortRawHitIgnoringDurationInTicks,
    kConsecutiveHitSeperationDurationInTicks,
    kGSFilterWindow,
    kGSFilterDegree,
    kCoincidenceWindowLowerLim,
    kCoincidenceWindowUpperLim
};

template <class T>
struct hit_t {
    T TStartInNanoSec;
    T TPeakInNanoSec;
    T IntegratedChargeInADCTimeTicks;
    T IntegratedChargeInADCNanoSec;
    T AmplitudeInMiliVolt;
    T AmplitudeInADC;
    T RiseTimeInNanoSec;
    T FallTimeInNanoSec;
    T CoincidenceWindowLowerLim;
    T CoincidenceWindowUpperLim;
    bool IsContained;
};

template <class T>
class CFDHitFinder {
public:
    CFDHitFinder() {
        _HitCollection.clear();
        _RawHitLogicNanosec.clear();
        _PedestalInADC = 0;
        _NoiseSigmaInADC = 0;
        _CFDThresholdInADC = 0;
        _WaveformADCNanosec.clear();
        _NonfilterWaveformADCNanosec.clear();
        _CFDParamSet.clear();
    }
    virtual ~CFDHitFinder() = default;

    void SetWaveform(std::vector<uint16_t>& waveform);
    void SetParam(CFDParams param, T value);

    inline const std::map<T, hit_t<T> >& GetHitCollection() const { return _HitCollection; }
    inline const T& GetPedestal() const { return _PedestalInADC; }
    inline const T& GetNoiseSigma() const { return _NoiseSigmaInADC; }
    inline const T& GetCFDThreshold() const { return _CFDThresholdInADC; }

    virtual void Go();
private:
    std::vector<T> _WaveformADCNanosec;
    std::vector<T> _NonfilterWaveformADCNanosec;
    std::vector<bool> _RawHitLogicNanosec;
    T _PedestalInADC;
    T _NoiseSigmaInADC;
    T _CFDThresholdInADC;

    std::map<CFDParams, T> _CFDParamSet;
    std::map<T, hit_t<T> > _HitCollection;

private:
    virtual void FindPedestal();
    virtual void FindRawHitLogic();
    virtual void FindCFDHits();
    virtual bool BackwardFindingOfHitStart(size_t hitPeakTimeIndex, T hitPeakValue, T& hitStartTimeIndex, T& hitRiseTimeInIndex);
    virtual bool ForwardFindingOfHitFallTime(size_t hitPeakTimeIndex, T& hitFallTimeInIndex);
    virtual T    IntegrateWaveformInADC(size_t hitStartTimeIndex, size_t hitEndTimeIndex);

    virtual inline void Reset() {
        _HitCollection.clear();
        _RawHitLogicNanosec.clear();
        _PedestalInADC = 0;
        _NoiseSigmaInADC = 0;
        _CFDThresholdInADC = 0;
        _WaveformADCNanosec.clear();
        _NonfilterWaveformADCNanosec.clear();
    }
};

template class hit_t<double>;

template class CFDHitFinder<double>;
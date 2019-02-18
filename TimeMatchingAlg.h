#pragma once

#include <vector>
#include <map>
#include <utility>

#include <HitFinderAlg.h>

template<class T>
class TimeMatchingAlg {
public:
    TimeMatchingAlg();
    virtual ~TimeMatchingAlg();

    void  SetCoincidenceTimeWindow(T coincidenceWindowLowerLim, T coincidenceWindowUpperLim, T timeSamplingInterval);
    virtual std::vector<std::pair<T, T>> GetAllTimeOfFlight(std::map<T, hit_t<T> > allCFDHitsInATriggerUpstream, std::map<T, hit_t<T> > allCFDHitsInATriggerDownstream);
    virtual std::vector<std::pair<T, T>> GetAllTimeOfFlightExperimental(std::map<T, hit_t<T> > allCFDHitsInATriggerUpstream, std::map<T, hit_t<T> > allCFDHitsInATriggerDownstream, T meanUS, T meanDS, T stdDevUS, T stdDevDS);

private:
    T _CoincidenceWindowLowerLim;
    T _CoincidenceWindowUpperLim;
    T _TimeSamplingInterval;
    std::vector<std::pair<T, T> > _AllTOFs;

private:
    std::vector<T> CheckInCoincidenceWindow(T hitStartUpstream, std::map<T, hit_t<T> > htStartDownstreamVector);
};

template class TimeMatchingAlg<double>;
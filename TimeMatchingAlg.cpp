#include <TimeMatchingAlg.h>

template<class T>
TimeMatchingAlg<T>::TimeMatchingAlg() : _CoincidenceWindowLowerLim(0), _CoincidenceWindowUpperLim(0), _TimeSamplingInterval(0)
{}

template<class T>
TimeMatchingAlg<T>::~TimeMatchingAlg() {}

template<class T>
std::vector<size_t> TimeMatchingAlg<T>::CheckInCoincidenceWindow(size_t hitStartUpstream, std::map<size_t, hit_t<T> > htStartDownstreamVector) {
    std::vector<size_t> timeDiffToThisUpstreamHit; timeDiffToThisUpstreamHit.clear();
    for (typename std::map<size_t, hit_t<T> >::iterator itrDS = htStartDownstreamVector.begin(); itrDS != htStartDownstreamVector.end(); itrDS++) {
        T timeDiff = (T)((*itrDS).first - hitStartUpstream) * _TimeSamplingInterval;
        if ((timeDiff - _CoincidenceWindowLowerLim >= 0) && (timeDiff - _CoincidenceWindowUpperLim <= 0)) {
            timeDiffToThisUpstreamHit.push_back((*itrDS).first);
        }
    }

    return timeDiffToThisUpstreamHit;
}


template<class T>
void TimeMatchingAlg<T>::SetCoincidenceTimeWindow(T coincidenceWindowLowerLim, T coincidenceWindowUpperLim, T timeSamplingInterval) {
    _CoincidenceWindowLowerLim = coincidenceWindowLowerLim;
    _CoincidenceWindowUpperLim = coincidenceWindowUpperLim;
    _TimeSamplingInterval = timeSamplingInterval;
}

template<class T>
std::vector<std::pair<T, T>> TimeMatchingAlg<T>::GetAllTimeOfFlight(std::map<size_t, hit_t<T> > allCFDHitsInATriggerUpstream, std::map<size_t, hit_t<T> > allCFDHitsInATriggerDownstream) {
    _AllTOFs.clear();

    for (typename std::map<size_t, hit_t<T> >::iterator itrUS = allCFDHitsInATriggerUpstream.begin(); itrUS != allCFDHitsInATriggerUpstream.end(); itrUS++) {
        std::vector<size_t> timeDiffToThisUpstreamHit = CheckInCoincidenceWindow(((*itrUS).first), allCFDHitsInATriggerDownstream);
        for (std::vector<size_t>::iterator itrCoincidence = timeDiffToThisUpstreamHit.begin(); itrCoincidence != timeDiffToThisUpstreamHit.end(); itrCoincidence++) {
            T chargeUS = itrUS->second.IntegratedChargeInADCNanoSec;
            T chargeDS = allCFDHitsInATriggerDownstream[*itrCoincidence].IntegratedChargeInADCNanoSec;
            std::pair<T, T> tmpPair = std::make_pair(((*itrCoincidence) - (*itrUS).first) * _TimeSamplingInterval,
                    (chargeUS - chargeDS));
            _AllTOFs.push_back(tmpPair);
        }
    }

    return _AllTOFs;
}

template<class T>
std::vector<std::pair<T, T>> TimeMatchingAlg<T>::GetAllTimeOfFlightExperimental(std::map<size_t, hit_t<T> > allCFDHitsInATriggerUpstream, std::map<size_t, hit_t<T> > allCFDHitsInATriggerDownstream, T meanUS, T meanDS, T stdDevUS, T stdDevDS) {

    _AllTOFs.clear();

    for (typename std::map<size_t, hit_t<T> >::iterator itrUS = allCFDHitsInATriggerUpstream.begin(); itrUS != allCFDHitsInATriggerUpstream.end(); itrUS++) {
        std::vector<size_t> timeDiffToThisUpstreamHit = CheckInCoincidenceWindow(((*itrUS).first), allCFDHitsInATriggerDownstream);
        for (std::vector<size_t>::iterator itrCoincidence = timeDiffToThisUpstreamHit.begin(); itrCoincidence != timeDiffToThisUpstreamHit.end(); itrCoincidence++) {
            T sigmaChargeUS = (itrUS->second.IntegratedChargeInADCNanoSec - meanUS) / stdDevUS;
            T sigmaChargeDS = (allCFDHitsInATriggerDownstream[*itrCoincidence].IntegratedChargeInADCNanoSec - meanDS) / stdDevDS;
            std::pair<T, T> tmpPair = std::make_pair(((*itrCoincidence) - (*itrUS).first) * _TimeSamplingInterval,
                                                     (sigmaChargeUS - sigmaChargeDS));
            _AllTOFs.push_back(tmpPair);
        }
    }

    return _AllTOFs;
}

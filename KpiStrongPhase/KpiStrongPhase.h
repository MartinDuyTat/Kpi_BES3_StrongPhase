// Martin Duy Tat 12th February 2021
/**
 * KpiStrongPhase is a BOSS algorithm
 * It runs a double tag analysis on the signal mode \f$D\to K_S^0K^+K^-\f$
 */

#ifndef KPISTRONGPHASE
#define KPISTRONGPHASE

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
// STL
#include<string>

class KpiStrongPhase: public Algorithm {
  public:
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KpiStrongPhase(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KpiStrongPhase();
    /**
     * This function runs when algorithm is initialized
     */
    StatusCode initialize();
    /**
     * Execution of the algorithm
     */
    StatusCode execute();
    /**
     * This function runs when algorithm is finalized
     */
    StatusCode finalize();
  private:
    /**
     * \f$K_SKK\f$ vs \f$K\pi\f$ tag mode
     */
    Algorithm *m_KSKKVersusKpiTag;
    /**
     * \f$K_SKK\f$ vs \f$K\pi\pi^0\f$ tag mode
     */
    Algorithm *m_KSKKVersusKpipi0Tag;
    /**
     * \f$K_SKK\f$ vs \f$K\pi\pi\pi\f$ tag mode
     */
    Algorithm *m_KSKKVersusKpipipiTag;
    /**
     * \f$K_SKK\f$ vs \f$Ke\nu\f$ tag mode
     */
    Algorithm *m_KSKKVersusKeNuTag;
    /**
     * \f$K_LKK\f$ vs \f$K\pi\f$ tag mode
     */
    Algorithm *m_KLKKVersusKpiTag;
    /**
     * Turn on \f$K_SKK\f$ vs \f$K\pi\f$ tag mode
     */
    bool m_recKSKKVersusKpiTag;
    /**
     * Turn on \f$K_SKK\f$ vs \f$K\pi\pi^0\f$ tag mode
     */
    bool m_recKSKKVersusKpipi0Tag;
    /**
     * Turn on \f$K_SKK\f$ vs \f$K\pi\pi\pi\f$ tag mode
     */
    bool m_recKSKKVersusKpipipiTag;
    /**
     * Turn on \f$K_SKK\f$ vs \f$Ke\nu\f$ tag mode
     */
    bool m_recKSKKVersusKeNuTag;
    /**
     * Turn on \f$K_LKK\f$ vs \f$K\pi\f$ tag mode
     */
    bool m_recKLKKVersusKpiTag;
};

#endif

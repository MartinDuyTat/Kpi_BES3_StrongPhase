#include "GaudiKernel/DeclareFactoryEntries.h"
#include "KpiStrongPhase/KpiStrongPhase.h"
#include "KpiStrongPhase/KSKKVersusKpiDoubleTag.h"
#include "KpiStrongPhase/KSKKVersusKpipi0DoubleTag.h"
#include "KpiStrongPhase/KSKKVersusKpipipiDoubleTag.h"
#include "KpiStrongPhase/KSKKVersusKeNuDoubleTag.h"
#include "KpiStrongPhase/KLKKVersusKpiDoubleTag.h"

DECLARE_ALGORITHM_FACTORY(KpiStrongPhase)
DECLARE_ALGORITHM_FACTORY(KSKKVersusKpiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KSKKVersusKpipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KSKKVersusKpipipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KSKKVersusKeNuDoubleTag)
DECLARE_ALGORITHM_FACTORY(KLKKVersusKpiDoubleTag)

DECLARE_FACTORY_ENTRIES(KpiStrongPhase) {
  DECLARE_ALGORITHM(KpiStrongPhase)
}



#pragma once

#include "DecayModelQ2W.h"

namespace elSpectro {

inline void ConfigureGlobalQ2Bias(bool use = true,
                                  double center = 1.8,
                                  double width = 0.45,
                                  double strength = 2.5,
                                  double lowFloor = 0.005,
                                  double power = 1.0) {
  DecayModelQ2W::SetGlobalUseHighQ2Bias(use);
  if (use) {
    DecayModelQ2W::SetGlobalHighQ2Bias(center, width, strength);
    DecayModelQ2W::SetGlobalHighQ2BiasShape(power, lowFloor);
  }
}

} // namespace elSpectro

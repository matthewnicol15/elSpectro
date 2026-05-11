#pragma once

#include "DecayModelst.h"

namespace elSpectro {

class DecayModelst_NoPS final : public DecayModelst {
public:
  using DecayModelst::DecayModelst;

  double PhaseSpaceWeightSq(double W) override {
    // Preserve side-effects that update dynamic daughter masses.
    DecayModel::PhaseSpaceWeightSq(W);

    // Remove all phase-space weighting at this decay stage.
    return 1.0;
  }
};

} // namespace elSpectro

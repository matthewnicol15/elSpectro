#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBodyFlat.h"
#include "DecayModelst_NoPS.h"
#include "PhaseSpaceDecay.h"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 3;
TH1F hQ2("Q2", "Q2", 1000, -2, 12);
TH2F hQ2_gE("Q2_gE", "", 1000, -2, 12, 1000, -2, 12);
TH2F hQ2_hyperon("Q2_hyperon", "", 1000, -2, 12, 1000, -2, 12);
TH1D heE("eE", "eE", 1000, 0, 20);
TH1D heTh("eTh", "eTh", 1000, 0, 180);
TH1D hKp1Th("hKp1Th", "hKp1Th", 100, 0, 180);
TH1D hKp2Th("hKp2Th", "hKp2Th", 100, 0, 180);
TH1D hKpAcceptance("hKpAcceptance", "hKpAcceptance", 4, 0, 3);
TH1D hYTh("YTh", "YTh", 1000, 0, 180);
TH1F hW("W", "W", 1000, 2, 5);
TH1F hV("Vertex", "Vertex", 100, 0, 20);
TH1F ht("t", "t", 1000, -2, 12);
TH1F hgE("gE", "gE", 1000, -2, 12);
TH1F hgTh("gTh", "gTh", 1000, 0, 180);
TH1F hHyperon1Rec("hHyperon1Rec", "; Y* to #Lambda #pi^{-} Mass (GeV)", 1000, 1.0, 4);
TH2F hKp1Th_v_hyperonRec("hKp1Th_v_hyperonRec", "", 1000, 1.0, 4, 1000, 0, 180);
TH1F hPKp("PKp", "K+", 100, 0, 10);
TH1F hPHyperon("PHyperon", "#Lambda", 100, 0, 10);
TH1F hPLambda("PLambda", "P", 100, 0, 10);
TH1F hPPi0("PPi0", "#pi^{-}", 100, 0, 10);

void ed_KpLambda(double ebeamE, int nEvents, int fileno)
{
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  gRandom->SetSeed(0);

  // define e- beam, pdg =11 momentum = _beamP
  auto elBeam = initial(11, ebeamE);
  auto elbeam = elBeam->GetInteracting4Vector();

  // deuteron target at rest
  // Use CDBonn potential to give nucleon momentum
  auto qfTarget = initial(2212, 0, 1000010020, model(new NuclearBreakup(2212, 2112)), new QuasiFreeNucleon(QuasiFree::CDBonnMomentum())); // CDBonn potential

  // get 4-momentum of target nucleon
  auto nucleon = qfTarget->GetInteracting4Vector();

  const double Lambda_rest = 1.115683;
  const double Pi_rest = 0.134976;
  const double Kp_rest = 0.493677;

  // ed -> e' n K+ Y(9995, broad) -> lambda pi0
  double wmax = (*elbeam + *nucleon).M();
  std::cout << "wmax " << wmax << std::endl;

  // Y*(9995) -> lambda pi0
  const double thr1 = Lambda_rest + Pi_rest;
  double max1 = wmax - 1.0 * Kp_rest;
  if (max1 <= thr1)
    max1 = thr1 + 1E-3;

  auto dist1 = new DistTF1{TF1(Form("YStarBW_%d", (int)ebeamE), "TMath::BreitWigner(x,3.0,3.0)", thr1, max1)};
  mass_distribution(9995, new DistTF1{TF1("YStarMass", "TMath::BreitWigner(x,3.0,3.0)", thr1, max1)});
  auto hyperon1 = static_cast<DecayingParticle *>(particle(9995, model(new GenericModelst{dist1, {}, {3122, 111}})));

  // decay of pGamma* to K+ Y*
  auto pGammaStarDecay = static_cast<DecayModelst *>(model(new DecayModelst_NoPS{{hyperon1}, {321}}));

  //
  // create mesonex electroproduction of X + neutron
  // TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  double s_strength = 0.1;
  double t_strength = 0.9;
  double t_slope = 0.1;

  mesonex(elBeam, qfTarget, new DecayModelQ2W{0, pGammaStarDecay, new TwoBody_stu{s_strength, t_strength, t_slope, 0, 0}});

  // give limits to the detected electron
  auto production = dynamic_cast<ElectronScattering *>(generator().Reaction());
  auto q2wModel = dynamic_cast<DecayModelQ2W *>(production->Model());

  production->SetLimitTarRest_eThmin(4.5 * TMath::DegToRad());
  production->SetLimitTarRest_eThmax(40 * TMath::DegToRad());

  // get pointers to produced particles fror diagnostic histos
  auto Lambda = hyperon1->Model()->Product(0);
  auto Pi0 = hyperon1->Model()->Product(1);

  auto Kp = pGammaStarDecay->Product(1);
  auto electron = dynamic_cast<DecayModelQ2W *>(production->Model())->GetScatteredElectron();

  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------

  writer(new LundWriter{Form("/home/nics/work/York/elSpectro/scripts/outputs/ed_to_KpLambda_unbiased_%d_%d.dat", (int)ebeamE, fileno)});

  // initilase the generator, may take some time for making distribution tables
  initGenerator();

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------

  gBenchmark->Start("e");

  int Acceptance = 0;
  int percentage = nEvents / 100;
  int total = 0;

  for (int i = 0; i < nEvents; i++)
  {
    nextEvent();

    total++;
    // fill diagnostic histograms
    auto photon = *elbeam - electron->P4();
    double Q2 = -photon.M2();

    double W = (photon + *prbeam).M();
    double t = -1.0 * (photon - Kp->P4()).M2();


    if (i % percentage == 0)
      std::cout << "event number " << i << std::endl;

    Acceptance = 0;

    double Kp1Theta = Kp->P4().Theta() * TMath::RadToDeg();
    hKp1Th.Fill(Kp1Theta);

    if (Kp1Theta >= 5 && Kp1Theta <= 45)
    {
      Acceptance++;
    }

    hKpAcceptance.Fill(Acceptance);

    hQ2.Fill(Q2);
    hW.Fill(W);
    ht.Fill(t);
    hgE.Fill(photon.E());
    hQ2_gE.Fill(Q2, photon.E());

    auto elec = electron->P4();
    heTh.Fill(elec.Theta() * TMath::RadToDeg());
    heE.Fill(elec.E());

    auto hyperonRec = Lambda->P4() + Pi0->P4();
    hHyperon1Rec.Fill(hyperonRec.M());

    hQ2_hyperon.Fill(Q2, hyperonRec.M());

    hKp1Th_v_hyperonRec.Fill(hyperonRec.M(), Kp1Theta);

    hPKp.Fill(Kp->P4().P());
    hPHyperon.Fill(hyperon1->P4().P());
    hPLambda.Fill(Lambda->P4().P());
    hPPi0.Fill(Pi0->P4().P());

  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  std::cout << "rejection rate is " << 100. * (total - nEvents) / total << std::endl;

  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D *)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D *)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("/home/nics/work/York/elSpectro/scripts/outputs/ed_to_KpLambda_unbiased_%d_%d.root", (int)ebeamE, fileno), "recreate");
  // total ep cross section inputs
  if (hWdist)
    hWdist->Write();
  if (hGenWdist)
    hGenWdist->Write();

  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hHyperon1Rec.Write();
  hKp1Th_v_hyperonRec.Write();
  hPKp.Write();
  hPHyperon.Write();
  hPLambda.Write();
  hPPi0.Write();
  hKp1Th.Write();
  hKpAcceptance.Write();
  hQ2_gE.Write();
  hQ2_hyperon.Write();
  fout->Close();

  generator().Summary();
}

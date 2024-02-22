#include <fstream>

#include "../commons/TreeInfo.h"
using namespace std;

void ConvertLPairToLHE(const std::string& filename = "events.root",
                       const std::string& output_lhe = "gammagammamumu.lpair_inelel_pt15_7tev.lhe",
                       int num_events = -1,
                       int first_event = 0) {
  auto* f1 = TFile::Open(filename.c_str());

  lpair::TreeRun run;
  run.attach(f1);

  ofstream output(output_lhe);

  lpair::TreeEvent evt;
  evt.attach(f1, "h4444");
  if (evt.tree->GetEntriesFast() < 1)
    throw std::runtime_error("No event in the file!");
  evt.tree->GetEntry(0);

  output << "<LesHouchesEvents version=\"1.0\">" << endl;
  output << "<header>" << endl;
  output << "This file was created from the output of the LPAIR generator" << endl;
  output << "</header>" << endl;

  output << "<init>" << endl;
  output << Form(
      "%d %d %.10E %.10E 0 0 10042 10042 2 1\n", evt.pdg_id[0], evt.pdg_id[1], run.sqrt_s * 0.5, run.sqrt_s * 0.5);
  output << Form("%.10E %.10E % .10E %d\n", run.xsect, run.errxsect, 2.673112e-4, 0);
  output << "</init>" << endl;

  const unsigned long long max_event = (num_events > 0 ? num_events : evt.tree->GetEntriesFast());

  for (unsigned long long i = first_event; i < max_event; i++) {
    evt.tree->GetEntry(i);

    // JH - for filtering on e/mu/had tau decays
    int founde = 0;
    int foundmu = 0;
    int foundpik = 0;
    for (int k = 0; k < evt.np; k++) {
      if (fabs(evt.pdg_id[k]) == 13)
        foundmu++;
      if (fabs(evt.pdg_id[k]) == 11)
        founde++;
      if (fabs(evt.pdg_id[k]) == 211 || fabs(evt.pdg_id[k]) == 321)
        foundpik++;
    }

    /* For selecting e+mu tautau events */
    /*
	if(founde == 0 || foundmu == 0)
	continue;
      */

    /* For selecting lepton + 1-prong tautau events - probably only works for ElEl */
    /*
      if((founde > 0) && (foundpik != 1))
      	continue;
      if((foundmu > 0) && (foundpik != 1))
      	continue;
      if((foundmu == 0) && (founde == 0))
      	continue;
      if((foundmu > 0) && (founde > 0))
	continue;
      */

    output << "<event>" << endl;
    output << evt.np << "   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
    //	cout << "there are " << evt.np << " particles in this event\n";

    for (int j = 0; j < evt.np; j++) {
      TLorentzVector part;
      part.SetPtEtaPhiE(evt.pt[j], evt.eta[j], evt.phi[j], evt.E[j]);
      double pz;
      if (j < 2) {
        double mom_beam = sqrt(evt.E[j] * evt.E[j] - evt.m[j] * evt.m[j]);
        pz = (j == 0 ? +1 : -1) * mom_beam;
      } else
        //	Stupid trick to produce inelastic events in both directions!
        pz = i % 2 == 0 ? part.Pz() : -part.Pz();

      int status = 1, moth1 = evt.parent1[j], moth2 = evt.parent2[j];
      if (evt.role[j] == 1 || evt.role[j] == 2)
        status = -9;
      else if (evt.role[j] == 41) {
        status = -1;
        moth1 = 1;
      } else if (evt.role[j] == 42) {
        status = -1;
        moth1 = 2;
      } else if (evt.role[j] == 0) {
        status = -2;
      }

      output << Form("%d\t%d\t%d\t%d\t0\t0\t% .10E\t% .10E\t% .10E\t% .10E\t% .10E\t0.\t%g\n",
                     evt.pdg_id[j],
                     status,
                     moth1,
                     moth2,
                     part.Px(),
                     part.Py(),
                     pz,
                     evt.E[j],
                     evt.m[j],
                     1.);
      //	output << "P "<< i*N+j << " " << evt.pdg_id[j] << " " << px[j] << " " << py[j] << " " << pz[j] << " " << evt.E[j] << " 1 0 0 0 0" << endl;
    }

    output << "</event>" << endl;
  }
  output << "</LesHouchesEvents>" << endl;
  output.close();

  cout << "Converted " << max_event << " events" << endl;
}

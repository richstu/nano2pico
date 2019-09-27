#ifndef H_BTAG_WEIGHTER
#define H_BTAG_WEIGHTER

#include <string>

#include "TH3D.h"

#include "pico_tree.hpp"
#include "BTagEntry.hpp"
#include "BTagCalibration.hpp"
#include "BTagCalibrationReader.hpp"

class BTagWeighter{
public:

  explicit BTagWeighter(bool is_fast_sim = false, int year = 2016);

  double EventWeight(pico_tree &pico, BTagEntry::OperatingPoint op,
		     const std::string &bc_full_syst, const std::string &udsg_full_syst,
		     const std::string &bc_fast_syst, const std::string &udsg_fast_syst) const;

  double EventWeight(pico_tree &pico, BTagEntry::OperatingPoint op,
		     const std::string &bc_full_syst, const std::string &udsg_full_syst) const;

  double EventWeight(pico_tree &pico, const std::vector<BTagEntry::OperatingPoint> &ops,
		     const std::string &bc_full_syst, const std::string &udsg_full_syst) const;

  double EventWeight(pico_tree &pico, const std::vector<BTagEntry::OperatingPoint> &ops,
		     const std::string &bc_full_syst, const std::string &udsg_full_syst,
		     const std::string &bc_fast_syst, const std::string &udsg_fast_syst) const;

  double JetBTagWeight(pico_tree &pico, std::size_t ijet, BTagEntry::OperatingPoint op,
		       const std::string &bc_full_syst, const std::string &udsg_full_syst,
		       const std::string &bc_fast_syst, const std::string &udsg_fast_syst) const;

  double JetBTagWeight(pico_tree &pico, std::size_t ijet, BTagEntry::OperatingPoint op,
		       const std::string &bc_full_syst, const std::string &udsg_full_syst) const;

  double JetBTagWeight(pico_tree &pico, std::size_t ijet, const std::vector<BTagEntry::OperatingPoint> &ops,
		       const std::string &bc_full_syst, const std::string &udsg_full_syst) const;

  double JetBTagWeight(pico_tree &pico, std::size_t ijet, const std::vector<BTagEntry::OperatingPoint> &ops,
		       const std::string &bc_full_syst, const std::string &udsg_full_syst,
		       const std::string &bc_fast_syst, const std::string &udsg_fast_syst) const;

private:
  double GetMCTagEfficiency(int pdgId, float pT, float eta, BTagEntry::OperatingPoint op) const;

  static const std::vector<BTagEntry::OperatingPoint> op_pts_;
  static const std::vector<BTagEntry::JetFlavor> flavors_;

  std::unique_ptr<BTagCalibration> calib_deep_full_;
  std::unique_ptr<BTagCalibration> calib_deep_fast_;
  std::map<BTagEntry::OperatingPoint, std::unique_ptr<BTagCalibrationReader> > readers_deep_full_;
  std::map<BTagEntry::OperatingPoint, std::unique_ptr<BTagCalibrationReader> > readers_deep_fast_;
  std::vector<TH3D> btag_efficiencies_deep_;

  double deep_csv_loose_, deep_csv_medium_, deep_csv_tight_;

  bool is_fast_sim_;
};

#endif

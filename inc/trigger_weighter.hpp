#ifndef H_TRIGGER_WEIGHTER
#define H_TRIGGER_WEIGHTER

#include <string>
#include <utility>
#include <vector>

#include "correction.hpp"
#include "pico_tree.hpp"

enum class LeptonHLTStatus {fail_all=0, pass_lowerdilep=1, pass_upperdilep=2, pass_singlelep=3};

/*! \class TriggerWeighter
  
  \brief Utility to calculate trigger scale factors for the H->Zgamma analysis
 */
class TriggerWeighter{
public: 

  /*!\brief creates a TriggerWeighter object and loads in the JSON files needed
    for a particular year
  
    \param[in] year  year of data taking
   */
  TriggerWeighter(std::string year);

  /*!\brief returns MC-data scale factor for event in the format {value, systup, systdown}
     where the variations are obtained by taking value*syst(up|down)
 
  \param[in] pico pico n-tuple. Lepton and trigger branches must be filled
 */
  std::vector<float> GetSF(pico_tree &pico);

private:
  /*!\brief returns MC-data scale factor for event in the format {value, systup, systdown}
     where the variations are obtained by taking value*syst(up|down)
   
    \param[in] electron_pt list of electron pts
    \param[in] muon_pt list of muon pts
    \param[in] electron_eta list of electron etas
    \param[in] muon_eta list of muon etas
    \param[in] electron_phi list of electron phis
    \param[in] pass_singleel if event passes single electron trigger(s)
    \param[in] pass_singlemu if event passes single muon trigger(s)
    \param[in] pass_diel if event passes double electron trigger(s)
    \param[in] pass_dimu if event passes double muon trigger(s)
   */
  std::vector<float> GetSF(std::vector<float> electron_pt, 
      std::vector<float> muon_pt, std::vector<float> electron_eta, 
      std::vector<float> muon_eta, std::vector<float> electron_phi,
      bool pass_singleel, bool pass_singlemu, bool pass_diel, bool pass_dimu);

  /*!\brief returns probability (efficiency) for event to pass electron OR muon 
    triggers in the format {nominal value, up variation, down variation}
   
    \param[in] electron_pt list of electron pts
    \param[in] muon_pt list of muon pts
    \param[in] electron_eta list of electron etas
    \param[in] electron_phi list of electron phis
    \param[in] muon_eta list of muon etas
    \param[in] pass_singleel if event passes single electron trigger(s)
    \param[in] pass_singlemu if event passes single muon trigger(s)
    \param[in] pass_diel if event passes double electron trigger(s)
    \param[in] pass_dimu if event passes double muon trigger(s)
    \param[in] is_data sets whether data or MC probability is calculated
   */
  std::vector<float> GetTotalProbability(
      std::vector<float> electron_pt, std::vector<float> muon_pt, 
      std::vector<float> electron_eta, std::vector<float> muon_eta, 
      std::vector<float> electron_phi,
      bool pass_singleel, bool pass_singlemu, bool pass_diel, bool pass_dimu, 
      bool is_data);

  /*!\brief returns probability (efficiency) for event to pass single or dilepton
    triggers for a particular lepton flavor in the format 
    {nominal value, up variation, down variation}
   
    \param[in] lepton_pt list of lepton pts
    \param[in] lepton_eta list of lepton etas
    \param[in] lepton_phi list of lepton phis
    \param[in] pass_singlelep if event passes single lepton trigger(s)
    \param[in] pass_dilep if event passes single lepton trigger(s)
    \param[in] is_data sets whether data or MC probability is calculated
    \param[in] is_electron sets whether e or mu probability is calculated
   */
  std::vector<float> GetFlavorProbability(
      std::vector<float> lepton_pt, std::vector<float> lepton_eta, 
      std::vector<float> lepton_phi,
      bool pass_singlelep, bool pass_dilep, bool is_data, bool is_electron);

  /*!\brief returns probability (efficiency) for a lepton to pass a given trigger
    leg in the format {nominal value, up variation, down variation}
   
    \param[in] lepton_pt lepton pt
    \param[in] lepton_eta lepton eta
    \param[in] lepton_phi lepton phi
    \param[in] is_data sets whether data or MC probability is calculated
    \param[in] is_electron sets whether e or mu probability is calculated
    \param[in] trigger_leg sets which trigger is evaluated
   */
  std::vector<float> GetLeptonProbability(float lepton_pt, float lepton_eta, 
      float lepton_phi, bool is_data, bool is_electron, 
      LeptonHLTStatus trigger_leg);

  std::unique_ptr<correction::CorrectionSet> cs_ello_;
  std::unique_ptr<correction::CorrectionSet> cs_elup_;
  std::unique_ptr<correction::CorrectionSet> cs_elsi_;
  std::unique_ptr<correction::CorrectionSet> cs_ello_hole_;
  std::unique_ptr<correction::CorrectionSet> cs_elup_hole_;
  std::unique_ptr<correction::CorrectionSet> cs_elsi_hole_;
  std::unique_ptr<correction::CorrectionSet> cs_mulo_;
  std::unique_ptr<correction::CorrectionSet> cs_muup_;
  std::unique_ptr<correction::CorrectionSet> cs_musi_;
  correction::Correction::Ref map_diel_lower_dataeff_; 
  correction::Correction::Ref map_diel_lower_dataunc_;
  correction::Correction::Ref map_diel_lower_mceff_; 
  correction::Correction::Ref map_diel_lower_mcunc_; 
  correction::Correction::Ref map_diel_upper_dataeff_; 
  correction::Correction::Ref map_diel_upper_dataunc_; 
  correction::Correction::Ref map_diel_upper_mceff_; 
  correction::Correction::Ref map_diel_upper_mcunc_; 
  correction::Correction::Ref map_singleel_dataeff_; 
  correction::Correction::Ref map_singleel_dataunc_; 
  correction::Correction::Ref map_singleel_mceff_; 
  correction::Correction::Ref map_singleel_mcunc_; 
  correction::Correction::Ref map_hole_diel_lower_dataeff_; 
  correction::Correction::Ref map_hole_diel_lower_dataunc_;
  correction::Correction::Ref map_hole_diel_lower_mceff_; 
  correction::Correction::Ref map_hole_diel_lower_mcunc_; 
  correction::Correction::Ref map_hole_diel_upper_dataeff_; 
  correction::Correction::Ref map_hole_diel_upper_dataunc_; 
  correction::Correction::Ref map_hole_diel_upper_mceff_; 
  correction::Correction::Ref map_hole_diel_upper_mcunc_; 
  correction::Correction::Ref map_hole_singleel_dataeff_; 
  correction::Correction::Ref map_hole_singleel_dataunc_; 
  correction::Correction::Ref map_hole_singleel_mceff_; 
  correction::Correction::Ref map_hole_singleel_mcunc_; 
  correction::Correction::Ref map_dimu_lower_dataeff_; 
  correction::Correction::Ref map_dimu_lower_dataunc_; 
  correction::Correction::Ref map_dimu_lower_mceff_; 
  correction::Correction::Ref map_dimu_lower_mcunc_; 
  correction::Correction::Ref map_dimu_upper_dataeff_; 
  correction::Correction::Ref map_dimu_upper_dataunc_; 
  correction::Correction::Ref map_dimu_upper_mceff_; 
  correction::Correction::Ref map_dimu_upper_mcunc_; 
  correction::Correction::Ref map_singlemu_dataeff_; 
  correction::Correction::Ref map_singlemu_dataunc_; 
  correction::Correction::Ref map_singlemu_mceff_; 
  correction::Correction::Ref map_singlemu_mcunc_; 
  bool post_bpix_;
};

#endif

# Validate code
source set_env.sh
scons
./scripts/produce_unit_test_cross_sections.py --output_log unit_test_cross_section.log
./scripts/organize_missing_cross_section.py --log unit_test_cross_section.log
./scripts/validate_unit_test_cross_section.py --output_filename validate_unit_test_cross_section.log --golden_cross_section_log OLD_CODE/unit_test_cross_section.log --validate_cross_section_log NEW_CODE/unit_test_cross_section.log

./scripts/produce_unit_test_htozgamma_NanoAODv9.py --output_folder unit_test_htozgamma_nanoaodv9 --output_log unit_test_htozgamma_nanoaodv9.log
./scripts/validate_unit_test_picos.py --output_log_filename validate_unit_test_htozgamma_nanoaodv9.log --unit_test_log_filename unit_test_htozgamma_nanoaodv9.log --golden_base_folder OLD_CODE/unit_test_htozgamma_nanoaodv9 --validate_base_folder NEW_CODE/unit_test_htozgamma_nanoaodv9

# Tag code
git tag
git tag htozgamma_kingscanyon_v0
git push origin htozgamma_kingscanyon_v0

# Produce
./scripts/produce_htozgamma_picos.py -t htozgamma_kingscanyon_v0 -y 2016APV,2016,2017,2018 -n NanoAODv9 -b /net/cms11/cms11r0/pico --use_telegram
./scripts/produce_htozgamma_picos.py -t htozgamma_kingscanyon_v0 -y 2022,2022EE -n NanoAODv11 -b /net/cms11/cms11r0/pico --use_telegram
./scripts/produce_htozgamma_picos.py -t htozgamma_kingscanyon_v0 -y 2023 -n NanoAODv11p9 -b /net/cms11/cms11r0/pico --use_telegram

# Modify tag
git tag -d htozgamma_kingscanyon_v0
git push origin --delete htozgamma_kingscanyon_v0
git commit
git tag htozgamma_kingscanyon_v0
git push origin htozgamma_kingscanyon_v0

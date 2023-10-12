source set_env.sh
scons
./scripts/produce_unit_test_cross_sections.py --output_log unit_test_cross_section.log
./scripts/organize_missing_cross_section.py --log unit_test_cross_section.log
./scripts/validate_unit_test_cross_section.py --output_filename validate_unit_test_cross_section.log --golden_cross_section_log OLD_CODE/unit_test_cross_section.log --validate_cross_section_log NEW_CODE/unit_test_cross_section.log

./scripts/produce_unit_test_htozgamma_NanoAODv9.py --output_folder unit_test_htozgamma_nanoaodv9 --output_log unit_test_htozgamma_nanoaodv9.log
./scripts/validate_unit_test_picos.py --output_log_filename validate_unit_test_htozgamma_nanoaodv9.log --unit_test_log_filename unit_test_htozgamma_nanoaodv9.log --golden_base_folder OLD_CODE/unit_test_htozgamma_nanoaodv9 --validate_base_folder NEW_CODE/unit_test_htozgamma_nanoaodv9

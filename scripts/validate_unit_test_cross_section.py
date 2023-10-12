#!/usr/bin/env python3
import re
import argparse

# cross_sections[path] = (cross section (pb), year)
def load_cross_sections(cross_section_path):
  cross_sections = {}
  with open(cross_section_path) as cross_section_file:
    for line in cross_section_file:
      path, year, cross_section = re.findall("filepath: (.*) year: (\d*|\d*APV|\d*EE) cross-section: (.*) pb",line)[0]
      if path in cross_sections: print("[Error] "+path+" is already inside cross_sections")
      cross_sections[path] = [cross_section, year]
  return cross_sections

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='''\Compares cross section logs.''', formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-o','--output_filename', required=True, help='Contains list of NanoAODs that have different cross sections.')
  parser.add_argument('-g','--golden_cross_section_log', required=True, help='Path for golden cross_section log')
  parser.add_argument('-v','--validate_cross_section_log', required=True, help='Path for validate cross_section log')
  args = parser.parse_args()

  #golden_cross_section_path = "/homes/jbkim/analysis/nano2pico/unit_test_cross_section.log"
  #validate_cross_section_path = "/homes/jbkim/analysis/nano2pico.variables/unit_test_cross_section.log"
  golden_cross_section_path = args.golden_cross_section_log
  validate_cross_section_path = args.validate_cross_section_log
  
  # cross_sections[path] = (cross section (pb), year)
  golden_cross_sections = load_cross_sections(golden_cross_section_path)
  validate_cross_sections = load_cross_sections(validate_cross_section_path)

  output_log = open(args.output_filename,'w')

  # Compare cross_sections
  is_different = False
  for path in golden_cross_sections:
    golden_cross_section = golden_cross_sections[path][0]
    if path not in validate_cross_sections:
      print(path+' missing in validation')
      output_log.write('missing '+path+' in validation\n')
      continue
    validate_cross_section = validate_cross_sections[path][0]
    if golden_cross_section != validate_cross_section:
      if is_different == False: 
        print(golden_cross_section_path+" and "+validate_cross_section_path+" have different cross sections.")
        is_different = True
      print(path+" "+golden_cross_section+" vs "+validate_cross_section+" (pb)")
      output_log.write(path+" "+golden_cross_section+" vs "+validate_cross_section+" (pb)\n")

  output_log.close()

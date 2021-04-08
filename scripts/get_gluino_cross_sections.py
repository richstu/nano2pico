#!/bin/env python
import json
import pandas as pd

# Download gluino json file from 
# wget https://raw.githubusercontent.com/fuenfundachtzig/xsec/master/json/pp13_gluino_NNLO%2BNNLL.json 
# ./scripts/get_gluino_cross_sections.py path_to_pp13_gluino_NNLO+NNLL.json

def loadJsonFiles(jsonFilenames, crossSectionData):
  # load data
  for jsonFilename in jsonFilenames:
    # Get input data
    data = json.load(open(jsonFilename))
    df   = pd.DataFrame.from_dict(data["data"], orient = "index")
    # restore mass as column and sort
    df["mass_GeV"] = df.index.astype(int)
    df = df.sort_values("mass_GeV")
    df.reset_index(inplace = True, drop = True)
    # Set output
    tag = data['process_id']
    crossSectionData[tag] = [[], [], []]
    for iPoint in range(len(df.mass_GeV)):
      mass = df.mass_GeV[iPoint]
      crossSection = df.xsec_pb[iPoint]
      crossSectionUncertainty = df.unc_pb[iPoint]
      crossSectionData[tag][0].append(mass)
      crossSectionData[tag][1].append(crossSection)
      crossSectionData[tag][2].append(crossSectionUncertainty)
    # Set label
    label = 'no label'
    if tag == 'pp13_glsq': label = '$\\tilde g\\tilde q; m_{\\tilde g} = m_{\\tilde q (u,d,c,s)} \\ll m_{\\tilde q (b,t)}$'
    elif tag == 'pp13_sqsq': label = '$\\tilde q\\tilde q^*; m_{\\tilde q} = m_{\\tilde q (u,d,c,s,b)} \\ll m_{\\tilde g, \\tilde t}$'
    elif tag == 'pp13_glgl': label = '$\\tilde g\\tilde g; m_{\\tilde g} \\ll m_{\\tilde q}$'
    elif tag == 'pp13_hino': label = '$\\tilde\\chi\\tilde\\chi; m_{\\tilde {\\chi}^\\pm}=m_{\\tilde {\\chi}^0}$'
    elif tag == 'pp13_stsb': label = '$\\tilde t\\tilde t^*; m_{\\tilde t} \\ll m_{\\tilde g}, m_{\\tilde q (u,d,c,s,b)}$'
    crossSectionData[tag].append(label)

if __name__ == "__main__":

  # crossSectionData[tag] = [[mass [GeV]], [cross_section [pb]], [cross_section_uncertainty [pb]], label]
  crossSectionData = {}
  jsonFilenames = ['pp13_gluino_NNLO+NNLL.json',] # pp13_glgl: ~g~g; decoupled:~q
  loadJsonFiles(jsonFilenames, crossSectionData)
  print("// xsec is pb, xsec_unc is relative uncertainty")
  print('void gluinoCrossSection(int glu_mass, double &xsec, double &xsec_unc) {')
  for iPoint in xrange(len(crossSectionData['pp13_glgl'][0])):
    mass = crossSectionData['pp13_glgl'][0][iPoint]
    cross_section = crossSectionData['pp13_glgl'][1][iPoint]
    cross_section_uncertainty = crossSectionData['pp13_glgl'][2][iPoint]
    if iPoint==0:
      print('  if (glu_mass == '+str(mass)+') { xsec = '+str(cross_section)+'; xsec_unc = '+str(cross_section_uncertainty/cross_section)+'; return;}')
    else: 
      print('  else if (glu_mass == '+str(mass)+') { xsec = '+str(cross_section)+'; xsec_unc = '+str(cross_section_uncertainty/cross_section)+'; return;}')
  print('  else { xsec = 0; xsec_unc = 0; }')
  print('}')
  

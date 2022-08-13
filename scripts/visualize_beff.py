#!/usr/bin/env python
import sys
import ROOT
import array
import argparse

ROOT.gROOT.SetBatch(True)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Draws efficiency", formatter_class=argparse.RawTextHelpFormatter)

  parser.add_argument('-i','--input_filename', required=True, help='Root file')
  parser.add_argument('-n','--hist_name', help='Histogram filename')
  parser.add_argument('-o','--output_folder', required=True, help='Folder to hold output pdfs')
  parser.add_argument('-p','--pdf_filename', required=True, help='Output pdf base filename')
  
  args = parser.parse_args()

  # Input
  #eff_filename = input_folder+"/btagEfficiency_DeepCSV_2016.root"
  eff_filename = args.input_filename
  # Output
  output_folder = args.output_folder
  #eff_image_filename = output_folder+"/btagEfficiency_DeepCSV_2016"
  eff_image_filename = output_folder+"/"+args.pdf_filename

  print("With file: "+eff_filename)
  root_file = ROOT.TFile(eff_filename)

  if not args.hist_name:
    print("[Error] Select one of the histograms with -n option. keys in "+args.input_filename+" are as following")
    root_file.ls()
    sys.exit()

  #eff_hist_name = "btagEfficiency_DeepCSV_loose"
  eff_hist_name = args.hist_name
  eff_hist = root_file.Get(eff_hist_name)
  nbins_eta = eff_hist.GetNbinsX()
  nbins_pt = eff_hist.GetNbinsY()
  nbins_flavor = eff_hist.GetNbinsZ()

  # Print
  output_log_filename = eff_image_filename+".log"
  output_log = open(output_log_filename,'w')
  for iFlavor in range(nbins_flavor+2):
    for iEta in range(nbins_eta+2):
      for iPt in range(nbins_pt+2):
        if iFlavor == 0 or iFlavor == nbins_flavor+1: continue
        if iEta == 0 or iEta == nbins_eta+1: continue
        if iPt == 0 or iPt == nbins_pt+1: continue
        eff = eff_hist.GetBinContent(iEta, iPt, iFlavor)
        eta_low = eff_hist.GetXaxis().GetBinLowEdge(iEta)
        eta_high = eff_hist.GetXaxis().GetBinUpEdge(iEta)
        pt_low = eff_hist.GetYaxis().GetBinLowEdge(iPt)
        pt_high = eff_hist.GetYaxis().GetBinUpEdge(iPt)
        flavor_low = eff_hist.GetZaxis().GetBinLowEdge(iFlavor)
        flavor_high = eff_hist.GetZaxis().GetBinUpEdge(iFlavor)
        line = "eff iEta["+str(eta_low)+","+str(eta_high)+"] iPt["+str(pt_low)+","+str(pt_high)+"] iFlavor["+str(flavor_low)+","+str(flavor_high)+"]: "+str(eff)
        print(line)
        output_log.write(line+"\n")

  # Draw across flavor
  for iFlavor in range(nbins_flavor+2):
    if iFlavor == 0 or iFlavor == nbins_flavor+1: continue
    flavor_center = eff_hist.GetZaxis().GetBinCenter(iFlavor)
    graph_name = "eff_flavor_"+"{:.0f}".format(flavor_center)
    graph_x = array.array('d')
    graph_y = array.array('d')
    eta_lines = []
    pt_texts = []
    eta_text_content = "abs(eta)="
    # Unroll eta,pt graph
    for iPt in range(nbins_pt+2):
      for iEta in range(nbins_eta+2):
        iBin = iPt*(nbins_eta+2)+iEta
        #print("bin: "+str(iBin)+" eff: "+str(eff_hist.GetBinContent(iEta, iPt, iFlavor)))
        graph_x.append(iBin)
        graph_y.append(eff_hist.GetBinContent(iEta, iPt, iFlavor))
        # Make text for eta
        if iPt == 0:
          eta_low = eff_hist.GetXaxis().GetBinLowEdge(iEta)
          eta_high = eff_hist.GetXaxis().GetBinUpEdge(iEta)
          eta_text_content += "[{:.1f}".format(eta_low)+","+"{:.1f}".format(eta_high)+"] "
      # line for each eta
      eta_lines.append(ROOT.TLine(iBin+0.5, 0, iBin+0.5, 1))
      # text for each pt
      pt_low = eff_hist.GetYaxis().GetBinLowEdge(iPt)
      pt_high = eff_hist.GetYaxis().GetBinUpEdge(iPt)
      pt_texts.append(ROOT.TLatex(iBin-(nbins_eta+2)+nbins_eta/2, 1.0, "["+"{:.0f}".format(pt_low)+","+"{:.0f}".format(pt_high)+"]"))
    pt_texts.append(ROOT.TLatex(0, 1.03, "pt ranges"))
    flavor_graph = ROOT.TGraph(len(graph_x), graph_x, graph_y);
    flavor_graph.SetMarkerStyle(21)
    flavor_graph.SetTitle(graph_name)
    flavor_graph.SetMinimum(0)
    flavor_graph.SetMaximum(1)
    flavor_canvas = ROOT.TCanvas(graph_name, graph_name, 700, 500)
    flavor_graph.Draw("AP")
    for eta_line in eta_lines: eta_line.Draw()
    for pt_text in pt_texts: 
      pt_text.SetTextSize(0.023)
      pt_text.Draw()
    eta_text = ROOT.TLatex(0, -0.1, eta_text_content)
    eta_text.SetTextSize(0.04)
    eta_text.Draw()
    flavor_canvas.SaveAs(eff_image_filename+"_"+graph_name+".pdf")




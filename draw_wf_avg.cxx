#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace po = boost::program_options;

struct Entry {
  int channel_id;
  unsigned int copy_id;
  int num_events;
  int num_hg;
  int num_lg;
  bool is_nnvt;
  double theta;
  double phi;
  std::vector<double> waveform;
  std::vector<double> stddev;
};

struct SelectedEntry {
  Entry data;
  int   group_idx;  // index into targets[]
};

int main(int argc, char **argv) {
  std::string input_file  = "wf_avg.root";
  std::string output_file = "wf_avg_plots.pdf";
  int max_per_group       = 10;
  bool no_band            = false;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "print help")
    ("input,i", po::value<std::string>(&input_file)->default_value("wf_avg.root"),
     "input ROOT file")
    ("output,o", po::value<std::string>(&output_file)->default_value("wf_avg_plots.pdf"),
     "output PDF file")
    ("max-per-group,n", po::value<int>(&max_per_group)->default_value(10),
     "max channels per theta group")
    ("no-band", po::bool_switch(&no_band)->default_value(false),
     "disable ±1σ confidence band");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 0;
  }
  bool draw_band = !no_band;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile f(input_file.c_str());
  auto *t = (TTree *)f.Get("wf_avg");
  if (!t) {
    std::cerr << "Tree 'wf_avg' not found in " << input_file << '\n';
    return 1;
  }

  int    channel_id, num_events, num_hg, num_lg;
  unsigned int copy_id;
  bool is_nnvt;
  double theta, phi;
  std::vector<double> *waveform = nullptr;
  std::vector<double> *stddev   = nullptr;

  t->SetBranchAddress("channelId", &channel_id);
  t->SetBranchAddress("copyId",    &copy_id);
  t->SetBranchAddress("numEvents", &num_events);
  t->SetBranchAddress("numHG",     &num_hg);
  t->SetBranchAddress("numLG",     &num_lg);
  t->SetBranchAddress("theta",     &theta);
  t->SetBranchAddress("phi",       &phi);
  t->SetBranchAddress("isNNVT",    &is_nnvt);
  t->SetBranchAddress("waveform",  &waveform);
  t->SetBranchAddress("stddev",    &stddev);

  std::vector<Entry> entries;
  entries.reserve(t->GetEntries());

  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    entries.push_back({channel_id, copy_id, num_events, num_hg, num_lg,
                       is_nnvt, theta, phi, *waveform, *stddev});
  }
  f.Close();

  constexpr double pi          = M_PI;
  constexpr double targets[]   = {0.0, pi / 4.0, pi / 2.0,
                                  3.0 * pi / 4.0, pi};
  constexpr const char *tnames[] = {"0", "#pi/4", "#pi/2", "3#pi/4", "#pi"};
  constexpr int n_targets = sizeof(targets) / sizeof(targets[0]);
  constexpr double window = 0.1;

  std::vector<SelectedEntry> selected;

  for (int t_idx = 0; t_idx < n_targets; ++t_idx) {
    double target = targets[t_idx];

    std::vector<std::pair<double, size_t>> candidates;
    for (size_t i = 0; i < entries.size(); ++i) {
      double delta = std::fabs(entries[i].theta - target);
      if (delta < window)
        candidates.emplace_back(delta, i);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    int n_to_take = std::min(max_per_group, (int)candidates.size());

    std::vector<SelectedEntry> group;
    for (int k = 0; k < n_to_take; ++k)
      group.push_back({entries[candidates[k].second], t_idx});

    std::sort(group.begin(), group.end(),
              [](const SelectedEntry &a, const SelectedEntry &b) {
                return a.data.phi < b.data.phi;
              });

    selected.insert(selected.end(), group.begin(), group.end());
  }

  constexpr int nsamples = 1008;

  TCanvas c("c", "", 1000, 700);
  bool first_page = true;

  for (size_t i = 0; i < selected.size(); ++i) {
    const auto &e = selected[i].data;
    int group_tag = selected[i].group_idx;

    std::vector<double> x(nsamples);
    std::vector<double> wf(nsamples);
    std::vector<double> sd_upper(nsamples);
    std::vector<double> sd_lower(nsamples);

    double y_min = 1e99, y_max = -1e99;
    for (int j = 0; j < nsamples; ++j) {
      x[j] = j;
      wf[j] = e.waveform[j];
      if (wf[j] < y_min) y_min = wf[j];
      if (wf[j] > y_max) y_max = wf[j];
    }

    TGraph g_band;
    if (draw_band) {
      g_band.Set(nsamples * 2);
      for (int j = 0; j < nsamples; ++j) {
        double s = e.stddev[j];
        sd_lower[j] = wf[j] - s;
        sd_upper[j] = wf[j] + s;
        g_band.SetPoint(j, x[j], sd_lower[j]);
        g_band.SetPoint(nsamples * 2 - 1 - j, x[j], sd_upper[j]);
        if (sd_lower[j] < y_min) y_min = sd_lower[j];
        if (sd_upper[j] > y_max) y_max = sd_upper[j];
      }
      g_band.SetFillColorAlpha(kGray, 0.35);
      g_band.SetLineWidth(0);
    }

    double margin = (y_max - y_min) * 0.05;
    y_min -= margin;
    y_max += margin;

    TGraph g_wf(nsamples, x.data(), wf.data());
    g_wf.SetLineColor(kBlue);
    g_wf.SetLineWidth(2);

    c.cd();

    TH1F *hframe = c.DrawFrame(0, y_min, nsamples - 1, y_max);
    hframe->GetXaxis()->SetTitle("Time sample");
    hframe->GetYaxis()->SetTitle("ADC count");
    hframe->GetYaxis()->SetTitleOffset(1.4);

    if (draw_band)
      g_band.Draw("f same");
    g_wf.Draw("l same");

    TLatex latex;
    latex.SetNDC(true);
    latex.SetTextSize(0.035);

    latex.SetTextAlign(12);
    latex.DrawLatex(0.14, 0.92,
                    ("PMT copy #" + std::to_string(e.copy_id)
                     + " (" + (e.is_nnvt ? "N" : "H") + ")").c_str());

    char pos_str[128];
    if (e.num_lg == -1) {
      std::snprintf(pos_str, sizeof(pos_str),
                    "#theta = %.1f#circ   #phi = %.1f#circ   H=%d, LG removed",
                    e.theta * 180.0 / M_PI,
                    e.phi * 180.0 / M_PI,
                    e.num_hg);
    } else {
      std::snprintf(pos_str, sizeof(pos_str),
                    "#theta = %.1f#circ   #phi = %.1f#circ   H=%d L=%d",
                    e.theta * 180.0 / M_PI,
                    e.phi * 180.0 / M_PI,
                    e.num_hg, e.num_lg);
    }
    latex.SetTextAlign(32);
    latex.DrawLatex(0.86, 0.92, pos_str);

    latex.SetTextAlign(22);
    latex.SetTextSize(0.06);
    latex.SetTextColorAlpha(kGray, 0.15);
    latex.DrawLatex(0.5, 0.5, ("#theta group: " + std::string(tnames[group_tag])).c_str());

    if (first_page) {
      c.Print((output_file + "(").c_str(), "pdf");
      first_page = false;
    } else {
      c.Print(output_file.c_str(), "pdf");
    }

    if ((i + 1) % 20 == 0)
      std::cout << "  page " << (i + 1) << " / " << selected.size() << '\n';
  }

  if (!first_page)
    c.Print((output_file + ")").c_str(), "pdf");

  std::cout << "Done: " << selected.size() << " pages -> " << output_file << '\n';
  return 0;
}

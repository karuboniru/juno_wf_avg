#include <Event/CdWaveformEvt.h>
#include <Event/CdWaveformHeader.h>
#include <EvtNavigator/EvtNavHelper.h>
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
#include <Geometry/CdGeom.h>
#include <Geometry/IRecGeomSvc.hh>
#include <Geometry/PmtGeom.h>
#include <Identifier/IDService.h>
#include <Identifier/Identifier.h>
#include <RootWriter/RootWriter.h>
#include <SniperKernel/AlgBase.h>
#include <SniperKernel/AlgFactory.h>
#include <SniperKernel/SniperLog.h>
#include <SniperKernel/SniperPtr.h>

#include <TTree.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

static constexpr int kWfLength = 1008;

// ─── baseline: most frequent rounded ADC value in average waveform ───

static double find_baseline(const std::array<double, kWfLength> &avg) {
  std::map<long, int> counts;
  for (double v : avg)
    ++counts[std::lround(v)];
  long peak = 0;
  int max_count = 0;
  for (const auto &[val, cnt] : counts) {
    if (cnt > max_count) {
      max_count = cnt;
      peak = val;
    }
  }
  return static_cast<double>(peak);
}

// ─── per-gain band accumulation ───

struct GainAccum {
  std::array<double, kWfLength> sum{};
  std::array<double, kWfLength> sq_sum{};
  int count{0};
};

// ══════════════════════════════════════════════════════════════════

class WfAverage final : public AlgBase {
public:
  WfAverage(const std::string &name) : AlgBase(name) {
    declProp("EnableGainCorrection", m_enable_gain_corr = true);
    declProp("HighGainScale", m_hg_scale = 0.08);
    declProp("LowGainScale", m_lg_scale = 0.55);
  }

  WfAverage(WfAverage &&) = delete;
  WfAverage(const WfAverage &) = delete;
  WfAverage &operator=(WfAverage &&) = delete;
  WfAverage &operator=(const WfAverage &) = delete;

  ~WfAverage() final = default;

  // ═══ initialize ═══

  bool initialize() final {
    LogInfo << "WfAverage::initialize" << '\n';

    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
      LogError << "Failed to get RootWriter" << '\n';
      return false;
    }

    m_tree = svc->bookTree(*getParent(), "USER_OUTPUT/wf_avg", "wf_average");
    m_tree->Branch("channelId", &m_out_channel_id);
    m_tree->Branch("copyId",    &m_out_copy_id);
    m_tree->Branch("numEvents", &m_out_num_events);
    m_tree->Branch("numHG",     &m_out_num_hg);
    m_tree->Branch("numLG",     &m_out_num_lg);
    m_tree->Branch("theta",     &m_out_theta);
    m_tree->Branch("phi",       &m_out_phi);
    m_tree->Branch("waveform", &m_out_waveform);
    m_tree->Branch("stddev",   &m_out_stddev);

    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
    if (navBuf.invalid()) {
      LogError << "Failed to get NavBuffer" << '\n';
      return false;
    }
    m_buf = navBuf.data();

    SniperPtr<IRecGeomSvc> recGeomSvc(getParent(), "RecGeomSvc");
    if (recGeomSvc.invalid()) {
      LogError << "Failed to get RecGeomSvc" << '\n';
      return false;
    }
    m_cdGeom = recGeomSvc->getCdGeom();

    IDService *idServ = IDService::getIdServ();
    idServ->init();

    return true;
  }

  // ═══ execute — accumulate raw ADC, separated by gain range ═══

  bool execute() final {
    ++m_events_visited;

    auto *nav = m_buf->curEvt();
    if (!nav) return true;

    auto hdr = JM::getHeaderObject<JM::CdWaveformHeader>(nav);
    if (!(hdr && hdr->hasEvent())) {
      if (m_events_visited % 1000 == 0)
        LogInfo << "processed " << m_events_visited
                << " events (" << m_events_with_cd << " with CD waveform)"
                << '\n';
      return true;
    }

    auto *evt = dynamic_cast<JM::CdWaveformEvt *>(hdr->event());
    if (!evt) return true;

    ++m_events_with_cd;

    const auto &channels = evt->channelData();

    for (const auto &[pmtId, waveform] : channels) {
      const auto &adc = waveform->adc();
      if (adc.size() != kWfLength) {
        LogFatal << "Unexpected waveform length " << adc.size()
                 << " for PMT " << pmtId << " (expected " << kWfLength << ')'
                 << '\n';
        std::abort();
      }

      bool is_hg = waveform->isHighGain();
      auto &gd   = is_hg ? m_hg[pmtId] : m_lg[pmtId];

      for (int i = 0; i < kWfLength; ++i) {
        double v = static_cast<double>(adc[i]);
        gd.sum[i]    += v;
        gd.sq_sum[i] += v * v;
      }
      ++gd.count;
    }

    return true;
  }

  // ═══ finalize — post-hoc baselines → scale → pool → output ═══

  bool finalize() final {
    LogInfo << "WfAverage::finalize — events visited: " << m_events_visited
            << ", with CD waveform: " << m_events_with_cd << '\n';

    // collect union of all channels from both gain maps
    std::set<int> all_channels;
    for (const auto &[id, _] : m_hg) all_channels.insert(id);
    for (const auto &[id, _] : m_lg) all_channels.insert(id);

    int n_hg_only = 0, n_lg_only = 0, n_both = 0;

    const double ratio = m_hg_scale / m_lg_scale;

    for (int pmtId : all_channels) {
      auto it_hg = m_hg.find(pmtId);
      auto it_lg = m_lg.find(pmtId);

      const GainAccum *hg_acc = (it_hg != m_hg.end()) ? &it_hg->second : nullptr;
      const GainAccum *lg_acc = (it_lg != m_lg.end()) ? &it_lg->second : nullptr;

      // ── tally ──
      if (hg_acc && lg_acc) ++n_both;
      else if (hg_acc)      ++n_hg_only;
      else                  ++n_lg_only;

      // ── compute per-gain μ, σ ──
      std::array<double, kWfLength> mu_hg{}, sigma_hg{};
      std::array<double, kWfLength> mu_lg{}, sigma_lg{};
      int n_hg = 0, n_lg = 0;

      if (hg_acc) {
        n_hg = hg_acc->count;
        double n_inv = 1.0 / n_hg;
        for (int i = 0; i < kWfLength; ++i) {
          mu_hg[i]    = hg_acc->sum[i] * n_inv;
          double var  = hg_acc->sq_sum[i] * n_inv - mu_hg[i] * mu_hg[i];
          sigma_hg[i] = std::sqrt(std::max(0.0, var));
        }
      }
      if (lg_acc) {
        n_lg = lg_acc->count;
        double n_inv = 1.0 / n_lg;
        for (int i = 0; i < kWfLength; ++i) {
          mu_lg[i]    = lg_acc->sum[i] * n_inv;
          double var  = lg_acc->sq_sum[i] * n_inv - mu_lg[i] * mu_lg[i];
          sigma_lg[i] = std::sqrt(std::max(0.0, var));
        }
      }

      // ── high-gain correction (if enabled) ──
      if (m_enable_gain_corr && n_hg > 0) {
        double baseline = find_baseline(mu_hg);
        for (int i = 0; i < kWfLength; ++i) {
          mu_hg[i]    = (mu_hg[i] - baseline) * ratio + baseline;
          sigma_hg[i] = sigma_hg[i] * ratio;
        }
      }

      // ── weighted pool ──
      int n_total = n_hg + n_lg;
      std::array<double, kWfLength> mu{}, sigma{};

      if (n_hg == 0) {
        mu    = mu_lg;
        sigma = sigma_lg;
      } else if (n_lg == 0) {
        mu    = mu_hg;
        sigma = sigma_hg;
      } else {
        double inv_n = 1.0 / n_total;
        for (int i = 0; i < kWfLength; ++i) {
          mu[i] = (n_hg * mu_hg[i] + n_lg * mu_lg[i]) * inv_n;

          double between_var = n_hg * (mu_hg[i] - mu[i]) * (mu_hg[i] - mu[i])
                             + n_lg * (mu_lg[i] - mu[i]) * (mu_lg[i] - mu[i]);
          double within_var  = n_hg * sigma_hg[i] * sigma_hg[i]
                             + n_lg * sigma_lg[i] * sigma_lg[i];
          sigma[i] = std::sqrt((within_var + between_var) * inv_n);
        }
      }

      // ── fill output ──
      m_out_channel_id = pmtId;
      m_out_num_events = n_total;
      m_out_num_hg     = n_hg;
      m_out_num_lg     = n_lg;

      Identifier pid(pmtId);
      auto *idServ = IDService::getIdServ();
      m_out_copy_id = idServ->id2CopyNo(pid);
      auto *pmt = m_cdGeom->getPmt(pid);
      if (!pmt) {
        LogFatal << "Failed to get PMT geometry for channel " << pmtId << '\n';
        std::abort();
      }
      m_out_theta = pmt->getCenter().Theta();
      m_out_phi   = pmt->getCenter().Phi();

      m_out_waveform.assign(mu.begin(), mu.end());
      m_out_stddev.assign(sigma.begin(), sigma.end());

      m_tree->Fill();
    }

    LogInfo << "Channels: " << all_channels.size()
            << "  (HG-only: " << n_hg_only
            << ", LG-only: " << n_lg_only
            << ", both: "   << n_both << ")\n";

    return true;
  }

private:
  JM::NavBuffer *m_buf{};
  TTree         *m_tree{};
  CdGeom        *m_cdGeom{};

  std::map<int, GainAccum> m_hg;
  std::map<int, GainAccum> m_lg;

  int m_events_visited{0};
  int m_events_with_cd{0};

  bool   m_enable_gain_corr{};
  double m_hg_scale{};
  double m_lg_scale{};

  int                m_out_channel_id{};
  unsigned int       m_out_copy_id{};
  int                m_out_num_events{};
  int                m_out_num_hg{};
  int                m_out_num_lg{};
  double             m_out_theta{};
  double             m_out_phi{};
  std::vector<double> m_out_waveform;
  std::vector<double> m_out_stddev;
};

DECLARE_ALGORITHM(WfAverage);

#include <Event/CdWaveformEvt.h>
#include <Event/CdWaveformHeader.h>
#include <EvtNavigator/EvtNavHelper.h>
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
#include <Geometry/CdGeom.h>
#include <Geometry/IPMTParamSvc.h>
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
#include <numeric>
#include <vector>

static constexpr int kWfLength = 1008;

struct GainAccum {
  std::array<double, kWfLength> sum{};
  std::array<double, kWfLength> sq_sum{};
  int count{0};
  int num_hg{0};
  int num_lg{0};
};

// ══════════════════════════════════════════════════════════════════

class WfAverage final : public AlgBase {
public:
  WfAverage(const std::string &name) : AlgBase(name) {
    declProp("HighGainScale", m_hg_scale = 0.08);
    declProp("LowGainScale", m_lg_scale = 0.55);
    declProp("BaselineSampleCount", m_baseline_n = 100);
    declProp("IgnoreLowGain", m_ignore_lg = false);
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
    m_tree->Branch("isNNVT",    &m_out_is_nnvt);
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

    SniperPtr<IPMTParamSvc> pmtSvc(getParent(), "PMTParamSvc");
    if (pmtSvc.invalid()) {
      LogError << "Failed to get PMTParamSvc" << '\n';
      return false;
    }
    m_pmt_svc = pmtSvc.data();

    return true;
  }

  // ═══ execute — per-event baseline subtraction + inline gain scaling ═══

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
    const double ratio = m_hg_scale / m_lg_scale;

    for (const auto &[pmtId, waveform] : channels) {
      const auto &adc = waveform->adc();
      if (adc.size() != kWfLength) {
        LogFatal << "Unexpected waveform length " << adc.size()
                 << " for PMT " << pmtId << " (expected " << kWfLength << ')'
                 << '\n';
        std::abort();
      }

      bool is_hg = waveform->isHighGain();

      if (m_ignore_lg && !is_hg)
        continue;

      double baseline = std::accumulate(adc.begin(),
                                        adc.begin() + m_baseline_n, 0.0)
                      / m_baseline_n;

      double scale = is_hg ? ratio : 1.0;

      auto &gd = m_acc[pmtId];

      for (int i = 0; i < kWfLength; ++i) {
        double v = (static_cast<double>(adc[i]) - baseline) * scale;
        gd.sum[i]    += v;
        gd.sq_sum[i] += v * v;
      }
      ++gd.count;
      if (is_hg) ++gd.num_hg;
      else       ++gd.num_lg;
    }

    return true;
  }

  // ═══ finalize — compute μ, σ, fill output ═══

  bool finalize() final {
    LogInfo << "WfAverage::finalize — events visited: " << m_events_visited
            << ", with CD waveform: " << m_events_with_cd << '\n';

    for (auto &[pmtId, gd] : m_acc) {
      m_out_channel_id = pmtId;
      m_out_num_events = gd.count;
      m_out_num_hg     = gd.num_hg;

      if (m_ignore_lg)
        m_out_num_lg = -1;
      else
        m_out_num_lg = gd.num_lg;

      Identifier pid(pmtId);
      auto *idServ = IDService::getIdServ();
      m_out_copy_id = idServ->id2CopyNo(pid);
      m_out_is_nnvt = m_pmt_svc->isNNVT(m_out_copy_id);
      auto *pmt = m_cdGeom->getPmt(pid);
      if (!pmt) {
        LogFatal << "Failed to get PMT geometry for channel " << pmtId << '\n';
        std::abort();
      }
      m_out_theta = pmt->getCenter().Theta();
      m_out_phi   = pmt->getCenter().Phi();

      m_out_waveform.resize(kWfLength);
      m_out_stddev.resize(kWfLength);

      double n_inv = 1.0 / static_cast<double>(gd.count);
      for (int i = 0; i < kWfLength; ++i) {
        double mean = gd.sum[i] * n_inv;
        double var  = gd.sq_sum[i] * n_inv - mean * mean;
        m_out_waveform[i] = mean;
        m_out_stddev[i]   = std::sqrt(std::max(0.0, var));
      }

      m_tree->Fill();
    }

    LogInfo << "Channels: " << m_acc.size() << '\n';
    return true;
  }

private:
  JM::NavBuffer *m_buf{};
  TTree         *m_tree{};
  CdGeom        *m_cdGeom{};
  IPMTParamSvc  *m_pmt_svc{};

  std::map<int, GainAccum> m_acc;

  int m_events_visited{0};
  int m_events_with_cd{0};

  bool   m_ignore_lg{};
  double m_hg_scale{};
  double m_lg_scale{};
  int    m_baseline_n{};

  int                m_out_channel_id{};
  unsigned int       m_out_copy_id{};
  int                m_out_num_events{};
  int                m_out_num_hg{};
  int                m_out_num_lg{};
  double             m_out_theta{};
  double             m_out_phi{};
  bool               m_out_is_nnvt{};
  std::vector<double> m_out_waveform;
  std::vector<double> m_out_stddev;
};

DECLARE_ALGORITHM(WfAverage);

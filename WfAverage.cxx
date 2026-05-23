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
#include <cstdint>
#include <cstdlib>
#include <map>
#include <vector>

static constexpr int kWfLength = 1008;

static double find_baseline(const std::vector<std::uint16_t> &adc) {
  std::map<std::uint16_t, int> counts;
  for (auto v : adc)
    ++counts[v];
  std::uint16_t peak = 0;
  int max_count = 0;
  for (const auto &[val, cnt] : counts) {
    if (cnt > max_count) {
      max_count = cnt;
      peak = val;
    }
  }
  return static_cast<double>(peak);
}

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
    m_tree->Branch("theta",      &m_out_theta);
    m_tree->Branch("phi",        &m_out_phi);
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

  bool execute() final {
    auto *nav = m_buf->curEvt();
    if (!nav) return true;

    auto hdr = JM::getHeaderObject<JM::CdWaveformHeader>(nav);
    if (!hdr || !hdr->hasEvent()) return true;

    auto *evt = dynamic_cast<JM::CdWaveformEvt *>(hdr->event());
    if (!evt) return true;

    const auto &channels = evt->channelData();

    for (const auto &[pmtId, waveform] : channels) {
      const auto &adc = waveform->adc();
      if (adc.size() != kWfLength) {
        LogFatal << "Unexpected waveform length " << adc.size()
                 << " for PMT " << pmtId << " (expected " << kWfLength << ')'
                 << '\n';
        std::abort();
      }

      std::vector<double> adc_d(kWfLength);
      if (m_enable_gain_corr && waveform->isHighGain()) {
        double baseline = find_baseline(adc);
        double ratio    = m_hg_scale / m_lg_scale;
        for (int i = 0; i < kWfLength; ++i)
          adc_d[i] = (static_cast<double>(adc[i]) - baseline) * ratio + baseline;
      } else {
        for (int i = 0; i < kWfLength; ++i)
          adc_d[i] = static_cast<double>(adc[i]);
      }

      auto &sum = m_sum[pmtId];
      auto &sq  = m_sq_sum[pmtId];

      for (int i = 0; i < kWfLength; ++i) {
        sum[i] += adc_d[i];
        sq[i]  += adc_d[i] * adc_d[i];
      }
      ++m_count[pmtId];
    }

    ++m_events_processed;
    return true;
  }

  bool finalize() final {
    LogInfo << "WfAverage::finalize — events with CD waveform: "
            << m_events_processed << '\n';
    LogInfo << "Channels with waveform data: " << m_count.size() << '\n';

    for (const auto &[pmtId, n] : m_count) {
      const auto &sum = m_sum[pmtId];
      const auto &sq  = m_sq_sum[pmtId];

      m_out_channel_id = pmtId;
      m_out_num_events = n;

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

      m_out_waveform.resize(kWfLength);
      m_out_stddev.resize(kWfLength);

      double n_inv = 1.0 / static_cast<double>(n);
      for (int i = 0; i < kWfLength; ++i) {
        double mean     = sum[i] * n_inv;
        double variance = sq[i] * n_inv - mean * mean;
        m_out_waveform[i] = mean;
        m_out_stddev[i]   = std::sqrt(std::max(0.0, variance));
      }

      m_tree->Fill();
    }

    return true;
  }

private:
  JM::NavBuffer *m_buf{};
  TTree         *m_tree{};
  CdGeom        *m_cdGeom{};

  std::map<int, std::array<double, kWfLength>> m_sum;
  std::map<int, std::array<double, kWfLength>> m_sq_sum;
  std::map<int, int>                            m_count;

  int m_events_processed{0};

  bool   m_enable_gain_corr{};
  double m_hg_scale{};
  double m_lg_scale{};

  int                m_out_channel_id{};
  unsigned int       m_out_copy_id{};
  int                m_out_num_events{};
  double             m_out_theta{};
  double             m_out_phi{};
  std::vector<double> m_out_waveform;
  std::vector<double> m_out_stddev;
};

DECLARE_ALGORITHM(WfAverage);

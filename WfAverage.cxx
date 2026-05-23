#include <Event/CdWaveformEvt.h>
#include <Event/CdWaveformHeader.h>
#include <EvtNavigator/EvtNavHelper.h>
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
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
#include <vector>

static constexpr int kWfLength = 1008;

class WfAverage final : public AlgBase {
public:
  WfAverage(const std::string &name) : AlgBase(name) {}

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
    m_tree->Branch("numEvents", &m_out_num_events);
    m_tree->Branch("waveform", &m_out_waveform);
    m_tree->Branch("stddev", &m_out_stddev);

    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
    if (navBuf.invalid()) {
      LogError << "Failed to get NavBuffer" << '\n';
      return false;
    }
    m_buf = navBuf.data();

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

      auto &sum = m_sum[pmtId];
      auto &sq  = m_sq_sum[pmtId];

      for (int i = 0; i < kWfLength; ++i) {
        double v = static_cast<double>(adc[i]);
        sum[i] += v;
        sq[i]  += v * v;
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

  std::map<int, std::array<double, kWfLength>> m_sum;
  std::map<int, std::array<double, kWfLength>> m_sq_sum;
  std::map<int, int>                            m_count;

  int m_events_processed{0};

  int                m_out_channel_id{};
  int                m_out_num_events{};
  std::vector<double> m_out_waveform;
  std::vector<double> m_out_stddev;
};

DECLARE_ALGORITHM(WfAverage);

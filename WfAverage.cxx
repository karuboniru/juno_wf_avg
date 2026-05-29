#include <Event/CdTriggerHeader.h>
#include <Event/CdWaveformEvt.h>
#include <Event/CdWaveformHeader.h>
#include <EvtNavigator/EvtNavHelper.h>
#include <EvtNavigator/EvtNavigator.h>
#include <EvtNavigator/NavBuffer.h>
#include <EDMPath/EDMPath.hh>
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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

static constexpr int kWfLength = 1008;

static inline double safe_sqrt(double x) {
  return std::sqrt(std::max(0.0, x));
}

struct GainAccum {
  std::array<double, kWfLength> sum{};
  std::array<double, kWfLength> sq_sum{};
  int count{0};
  int num_hg{0};
  int num_lg{0};
  std::array<double, 8> fadc_baseline_sum{};
  int fadc_baseline_count{0};
};

// ══════════════════════════════════════════════════════════════════

class WfAverage final : public AlgBase {
public:
  WfAverage(const std::string &name) : AlgBase(name) {
    declProp("HighGainScale", m_hg_scale = 0.08);
    declProp("LowGainScale", m_lg_scale = 0.55);
    declProp("BaselineSampleCount", m_baseline_n = 100);
    declProp("PerFadcBaseline", m_per_fadc_baseline = true);
    declProp("IgnoreLowGain", m_ignore_lg = false);
    declProp("TriggerTypeFilter", m_trig_filter = std::string{});
    declProp("TriggerInclusive", m_trigger_inclusive = false);
    declProp("TimeAlign", m_time_align = false);
    declProp("SkipOnMissingRef", m_skip_on_missing = true);
    declProp("MonitorChannel", m_monitor_channel = 43303);
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

    if (!m_trig_filter.empty()) {
      bool matched = false;
      if (auto trig_hdr = JM::getHeaderObject<JM::CdTriggerHeader>(nav))
        if (trig_hdr->hasEvent()) {
          const auto &types = trig_hdr->event()->triggerType();
          if (m_trigger_inclusive) {
            for (const auto &t : types)
              if (t == m_trig_filter) { matched = true; break; }
          } else {
            matched = (types.size() == 1 && types[0] == m_trig_filter);
          }
        }
      if (!matched) return true;
    }

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
    ++m_total_trigger_events;

    const auto &channels = evt->channelData();
    const double ratio = m_hg_scale / m_lg_scale;

    int delta_t = 0;

    if (m_time_align) {
      auto mon_hdr = JM::getHeaderObject<JM::CdWaveformHeader>(
          nav, JM::Elec::CdMonitorWaveform::Path);
      if (mon_hdr && mon_hdr->hasEvent()) {
        auto *mon_evt = dynamic_cast<JM::CdWaveformEvt *>(mon_hdr->event());
        if (mon_evt) {
          auto mon_it = mon_evt->channelData().find(m_monitor_channel);
          if (mon_it != mon_evt->channelData().end()) {
            const auto &mon_adc = mon_it->second->adc();
            if (mon_adc.size() == kWfLength) {
              auto min_it = std::min_element(mon_adc.begin(), mon_adc.end());
              int cur_peak_idx = static_cast<int>(std::distance(mon_adc.begin(), min_it));

              if (!m_ref_peak_time.has_value()) {
                m_ref_peak_time = cur_peak_idx;
              } else {
                delta_t = *m_ref_peak_time - cur_peak_idx;
                m_delta_sum += delta_t;
                m_delta_sq_sum += delta_t * delta_t;
                ++m_delta_count;
              }
            }
          } else {
            ++m_skipped_ref_missing;
            if (m_skip_on_missing) return true;
          }
        }
      }
    }

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

      double scale = is_hg ? ratio : 1.0;

      auto &gd = m_acc[pmtId];

      int start = 0, end = kWfLength, offset = 0;
      if (m_time_align && delta_t != 0) {
        start  = std::max(0, -delta_t);
        end    = std::min(kWfLength, kWfLength - delta_t);
        offset = delta_t;
      }

      if (m_per_fadc_baseline) {
        for (int i = 0; i < m_baseline_n; ++i)
          gd.fadc_baseline_sum[i & 7] += adc[i] * scale;
        ++gd.fadc_baseline_count;

        for (int i = start; i < end; ++i) {
          double v = static_cast<double>(adc[i]) * scale;
          gd.sum[i + offset]    += v;
          gd.sq_sum[i + offset] += v * v;
        }
      } else {
        double baseline = std::accumulate(adc.begin(),
                                          adc.begin() + m_baseline_n, 0.0)
                        / m_baseline_n;
        for (int i = start; i < end; ++i) {
          double v = (static_cast<double>(adc[i]) - baseline) * scale;
          gd.sum[i + offset]    += v;
          gd.sq_sum[i + offset] += v * v;
        }
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
            << ", with CD waveform: " << m_events_with_cd
            << ", trigger events: " << m_total_trigger_events;
    if (m_time_align) {
      double dt_mean = m_delta_count > 0
          ? static_cast<double>(m_delta_sum) / m_delta_count : 0.0;
      double dt_std  = m_delta_count > 1
          ? safe_sqrt(
              static_cast<double>(m_delta_sq_sum) / m_delta_count
              - dt_mean * dt_mean)
          : 0.0;
      LogInfo << " [time-align: monitor_chan=" << m_monitor_channel
              << " ref_peak=" << m_ref_peak_time.value_or(-1)
              << " skipped_missing=" << m_skipped_ref_missing
              << " jitter (μ ± σ): " << dt_mean << " ± " << dt_std << " bins"
              << " (N=" << m_delta_count << ")]";
    }
    LogInfo << '\n';

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
      if (auto *pmt = m_cdGeom->getPmt(pid); pmt) {
        m_out_theta = pmt->getCenter().Theta();
        m_out_phi   = pmt->getCenter().Phi();
      } else {
        LogFatal << "Failed to get PMT geometry for channel " << pmtId << '\n';
        std::abort();
      }

      m_out_waveform.resize(kWfLength);
      m_out_stddev.resize(kWfLength);

      double n_inv = 1.0 / static_cast<double>(m_total_trigger_events);
      for (int i = 0; i < kWfLength; ++i) {
        double mean = gd.sum[i] * n_inv;
        double var  = gd.sq_sum[i] * n_inv - mean * mean;
        m_out_waveform[i] = mean;
        m_out_stddev[i]   = safe_sqrt(var);
      }

      if (m_per_fadc_baseline && gd.fadc_baseline_count > 0) {
        double inv = 1.0 / static_cast<double>(gd.fadc_baseline_count);
        for (int i = 0; i < kWfLength; ++i)
          m_out_waveform[i] -= gd.fadc_baseline_sum[i & 7] * inv;
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
  int m_total_trigger_events{0};

  bool        m_ignore_lg{};
  bool        m_trigger_inclusive{};
  double      m_hg_scale{};
  double m_lg_scale{};
  int    m_baseline_n{};
  bool   m_per_fadc_baseline{};
  std::string m_trig_filter{};

  bool m_time_align{};
  bool m_skip_on_missing{};
  std::optional<int> m_ref_peak_time{};
  int  m_monitor_channel{43303};
  int  m_skipped_ref_missing{0};
  long m_delta_sum{0};
  long m_delta_sq_sum{0};
  int  m_delta_count{0};

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

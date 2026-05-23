# WfAverage — CD Waveform Event-Averaging for JUNO

对 JUNO 探测器 rtraw 原始数据中 `/Event/CdWaveform` 的波形进行 event-by-event
平均，输出每个 PMT channel 的平均波形及标准差，并支持 multi-gain 修正。

## 文件结构

```
wf_average/
├── CMakeLists.txt       # CMake 构建：MODULE (WfAverage) + EXECUTABLE (draw_wf_avg)
├── WfAverage.cxx        # SNiPER AlgBase 算法 — 读 rtraw → 累加 → 输出
├── run.py               # Python 驱动脚本（多个输入文件 / --input-list）
├── draw_wf_avg.cxx      # ROOT 独立画图程序 — 读 wf_avg.root → 多页 PDF
└── README.md
```

## 依赖

- JUNO Offline Software (junosw), 版本 J26.1.x 或更高
- C++ 23 (g++-15)
- ROOT 6.30+
- SNiPER v2

## 构建

```bash
source ~/.junorc
cmake -B build -G Ninja
cmake --build build
```

产物：
- `build/lib/libWfAverage.so` — SNiPER 算法模块
- `build/bin/draw_wf_avg.exe` — 画图程序

---

## 1. 波形平均（WfAverage）

### 运行

```bash
python run.py --input file1.rtraw file2.rtraw ...  [--evtmax N] [--user-output wf_avg.root]

# 或通过文本文件列出输入路径（每行一个）：
python run.py --input-list files.txt
```

### 算法属性（可在 Python 中覆盖）

| 属性 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `EnableGainCorrection` | bool | `true` | 是否启用 high-gain 修正 |
| `HighGainScale` | double | `0.08` | HG scale 因子 |
| `LowGainScale` | double | `0.55` | LG scale 因子 |

其中 $r = \frac{\text{HighGainScale}}{\text{LowGainScale}} \approx 0.145$ 为缩放比。

在 Python 中覆盖示例：
```python
alg = task.createAlg("WfAverage")
alg.property("EnableGainCorrection").set(False)
alg.property("HighGainScale").set(0.07)
```

### 输出 ROOT 文件

TTree `USER_OUTPUT/wf_avg`（内部名 `wf_average`），每行一个 channel：

| Branch | 类型 | 说明 |
|---|---|---|
| `channelId` | int | JUNO 原始 PMT 标识符 |
| `copyId` | unsigned int | 人类友好的复制编号（0-17611） |
| `numEvents` | int | 该 channel 出现的事件总数 (= numHG + numLG) |
| `numHG` | int | High-gain 事件数 |
| `numLG` | int | Low-gain 事件数 |
| `theta` | double | PMT 天顶角 [rad] |
| `phi` | double | PMT 方位角 [rad] |
| `waveform` | vector\<double\> | 平均波形（1008 个时间 bin） |
| `stddev` | vector\<double\> | 标准差（1008 个时间 bin） |

### 算法流程

#### execute() — 累加阶段（不做任何修正）

```
对每个 rtraw event:
  跳过没有 CdWaveformHeader 的事件
  获取 CdWaveformEvt → 遍历所有 channel:
    若 isHighGain() → 累加原始 ADC 到 m_hg[pmtId]
    若 isLowGain()  → 累加原始 ADC 到 m_lg[pmtId]

每 1000 个 event 打印进度。
```

#### finalize() — 后处理阶段

对 `m_hg ∪ m_lg` 中的每个 channel：

1. **计算各 gain band 的统计量**：从累加的 $\Sigma\mathrm{adc}$ 和 $\Sigma\mathrm{adc}^2$ 计算各 gain 的逐 bin 均值和标准差：

$$\mu_{\mathrm{HG}}[i] = \frac{\Sigma\mathrm{adc}_{\mathrm{HG}}[i]}{N_{\mathrm{HG}}}, \qquad
  \sigma_{\mathrm{HG}}[i] = \sqrt{\max\!\left(0,\;
    \frac{\Sigma\mathrm{adc}^2_{\mathrm{HG}}[i]}{N_{\mathrm{HG}}} - \mu_{\mathrm{HG}}[i]^2
  \right)}$$

   （LG 同理）。

2. **High-gain 修正**（仅当 `EnableGainCorrection` 且该 channel 同时有 HG 和 LG 事件时执行）：

$$\textit{baseline} = \mathrm{mode}\!\big(\lfloor\mu_{\mathrm{HG}}[i] + 0.5\rfloor\big)_{i=0}^{1007}$$

$$\mu'_{\mathrm{HG}}[i] = \big(\mu_{\mathrm{HG}}[i] - \textit{baseline}\big) \cdot r + \textit{baseline}, \qquad
  \sigma'_{\mathrm{HG}}[i] = \sigma_{\mathrm{HG}}[i] \cdot r$$

   其中 $r = \dfrac{0.08}{0.55} \approx 0.145$（`HighGainScale / LowGainScale`）。

3. **加权池化**（若 HG 和 LG 同时存在）：

$$\begin{aligned}
N &= N_{\mathrm{HG}} + N_{\mathrm{LG}} \\
\mu[i] &= \frac{N_{\mathrm{HG}} \cdot \mu'_{\mathrm{HG}}[i] + N_{\mathrm{LG}} \cdot \mu_{\mathrm{LG}}[i]}{N} \\
\sigma^2[i] &= \frac{1}{N}\Big[
  \underbrace{N_{\mathrm{HG}} \cdot \sigma'^2_{\mathrm{HG}}[i] + N_{\mathrm{LG}} \cdot \sigma^2_{\mathrm{LG}}[i]}_{\text{组内方差}}
  + \underbrace{N_{\mathrm{HG}}\big(\mu'_{\mathrm{HG}}[i] - \mu[i]\big)^2}_{\text{组间方差}}
  + \underbrace{N_{\mathrm{LG}}\big(\mu_{\mathrm{LG}}[i] - \mu[i]\big)^2}_{\text{组间方差}}
\Big] \\
\sigma[i] &= \sqrt{\max(0,\;\sigma^2[i])}
\end{aligned}$$

   若只有单一 gain：直接输出该 gain 的原值，不应用修正。

### 关于 gain 修正

JUNO 电子学对 PMT 信号有两种增益范围：
- **High-gain (HR)**: 较高放大倍数（$\approx 1 / r = 0.55 / 0.08 \approx 6.875\times$），适合小信号
- **Low-gain (LR)**: 较低放大倍数，适合大信号

不同 event 的同一 channel 可能触发不同 gain。直接混合会导致 scale 不一致。
基线扣除-缩放方法将 HG 波形相对于 pedestal (baseline) 缩放到与 LG 一致的
幅度尺度（乘以 $r = \frac{0.08}{0.55}$），使两者可安全地加权平均。

**修正仅在有 HG+LG 混合的 channel 上生效**。若 channel 只有一种 gain，无需
修正，原值直接输出。

---

## 2. 画图程序（draw_wf_avg）

### 运行

```bash
./build/bin/draw_wf_avg.exe [input.root] [output.pdf] [max_per_group]
```

参数：
- `input.root` — WfAverage 输出的 ROOT 文件（默认 `wf_avg.root`）
- `output.pdf` — 输出多页 PDF（默认 `wf_avg_plots.pdf`）
- `max_per_group` — 每个 θ 组选取的 channel 数（默认 10）

### 选取逻辑

在 5 个天顶角目标值附近各选最多 max_per_group 个 channel：

| θ 目标 | 组标签 |
|---|---|
| 0 | 0 |
| π/4 | π/4 |
| π/2 | π/2 |
| 3π/4 | 3π/4 |
| π | π |

选取窗口：±0.1 rad。在每个窗口内取最靠近目标的 `max_per_group` 个 channel，
按 φ 排序后在 PDF 中按 φ 顺序排列。

### 每页布局

- **左上角**: `PMT copy #XXXX`
- **右上角**: `θ = XX.X°  φ = XX.X°  H=X L=Y`
- **中央水印**: `θ group: π/4`（半透明）
- **主体**: 蓝色实线 = 平均波形，灰色半透明填充 = ±1σ 置信带
- **坐标轴**: X = 时间采样点 [0, 1007]，Y = ADC 计数

---

## 3. 典型工作流

```bash
# 1. 环境设置
source ~/.junorc

# 2. 构建
cmake -B build -G Ninja && cmake --build build

# 3. 运行平均算法
python run.py --input /path/to/run_*.rtraw --user-output wf_avg.root

# 4. 检查输出（ROOT 交互）
root -l wf_avg.root
root [1] USER_OUTPUT->cd("wf_avg")
root [2] wf_average->Scan("copyId:numEvents:theta:phi:numHG:numLG")

# 5. 生成画图
./build/bin/draw_wf_avg.exe wf_avg.root wf_avg_plots.pdf 10

# 6. 查看 PDF
evince wf_avg_plots.pdf
```

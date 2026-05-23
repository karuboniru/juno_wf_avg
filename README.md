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
| `HighGainScale` | double | `0.08` | HG scale 因子 |
| `LowGainScale` | double | `0.55` | LG scale 因子 |
| `BaselineSampleCount` | int | `100` | 用于基线估计的前 N 个采样点 |
| `IgnoreLowGain` | bool | `false` | 若 `true`，跳过所有 LG channel |

其中 $r = \frac{\text{HighGainScale}}{\text{LowGainScale}} \approx 0.145$ 为缩放比。

在 Python 中覆盖示例：
```python
alg = task.createAlg("WfAverage")
alg.property("BaselineSampleCount").set(200)
alg.property("IgnoreLowGain").set(True)
```

### 输出 ROOT 文件

TTree `USER_OUTPUT/wf_avg`（内部名 `wf_average`），每行一个 channel：

| Branch | 类型 | 说明 |
|---|---|---|
| `channelId` | int | JUNO 原始 PMT 标识符 |
| `copyId` | unsigned int | 人类友好的复制编号（0-17611） |
| `numEvents` | int | 该 channel 出现的事件总数 |
| `numHG` | int | High-gain 事件数 |
| `numLG` | int | Low-gain 事件数（`IgnoreLowGain=true` 时为 -1） |
| `theta` | double | PMT 天顶角 [rad] |
| `phi` | double | PMT 方位角 [rad] |
| `waveform` | vector\<double\> | 基线归零后的平均波形（1008 个时间 bin） |
| `stddev` | vector\<double\> | 标准差（1008 个时间 bin） |

> 波形值以 0 为基线：表示相对于 pedestal 的信号幅度。HG 已被
> scale 至与 LG 可比的幅度。

### 算法流程

#### execute() — 逐 event 基线减除 + 内联缩放

```
对每个 rtraw event:
  跳过没有 CdWaveformHeader 的事件
  获取 CdWaveformEvt → 遍历所有 channel:

    若 IgnoreLowGain && !isHighGain() → 跳过

    baseline = mean(adc[0 .. BaselineSampleCount - 1])   ← 前 N 个采样均值
    scale    = isHighGain() ? (HighGainScale / LowGainScale) : 1.0

    对每个采样点 i ∈ [0, 1007]:
      value = (adc[i] − baseline) × scale
      累加 value 和 value² 到 m_acc[pmtId]

  每 1000 个 event 打印进度。
```

#### finalize() — 输出统计量

对 `m_acc` 中每个 channel：

1. 计算逐 bin 均值和标准差：

$$\mu[i] = \frac{\Sigma v[i]}{N}, \qquad
  \sigma[i] = \sqrt{\max\!\left(0,\;
    \frac{\Sigma v^2[i]}{N} - \mu[i]^2
  \right)}$$

2. 若 `IgnoreLowGain` 启用，`numLG` 输出为 -1（哨兵值）；否则输出实际计数。
3. 填入输出 tree。

### 关于 gain 修正

JUNO 电子学对 PMT 信号有两种增益范围：
- **High-gain (HR)**: 较高放大倍数（$\approx 1 / r = 0.55 / 0.08 \approx 6.875\times$），适合小信号
- **Low-gain (LR)**: 较低放大倍数，适合大信号

**逐 event 修正策略**：对每个 channel 在每个 event 中：
1. 取前 `BaselineSampleCount` 个采样点的均值作为该 event 的基线（pedestal estimate）
2. 将波形减去该基线，得到信号幅度
3. 若为 HG channel，将幅度乘以 $r = 0.08/0.55$，缩放到与 LG 一致的尺度
4. LG channel 仅做基线减除（scale = 1.0），不做额外缩放

这样所有 event 的波形都以 0 为基线对齐，且 HG 与 LG 可直接混合累加。

`IgnoreLowGain` 模式：若该 channel 出现的是 LG 事件，直接跳过不累加。

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

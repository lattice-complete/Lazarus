# Lazarus ← lazer/labrador C 遷移地圖

本文檔追蹤 C 代碼庫(IBM [lazer](https://github.com/lazer-crypto/lazer) 及其
[labrador](https://github.com/lazer-crypto/labrador) 子模塊,後者基於
[lattice-dogs/labrador](https://github.com/lattice-dogs/labrador) 的算術)到
Lazarus Rust workspace 的遷移狀態與設計決策。

## 總體策略

可移植純 Rust 優先:先要正確、可測、跨平台(本次地基層已在 arm64 上驗證),
AVX-512 / NEON 的 SIMD 路徑留作後續 feature flag。C 版約 64K 行
(lazer src ~35K + labrador ~29K),按層遷移,證明系統層建立在本地基層之上。

## 地基層映射(本輪完成)

| C 來源 | Rust 位置 | 狀態 | 備註 |
|---|---|---|---|
| labrador `data32.c` 模數定義 | `algebra::zq` | ✅ | q = 2^32−99(素數);修復了原 `pow` 溢出,新增 `inv`/`neg`/`centered`/`from_i64` |
| lazer `poly.c`, labrador `poly.c`/`polz.c` | `algebra::ring::poly::Poly<D>`(`Rq64 = Poly<64>`) | ✅ | 度數進類型;NTT 乘法、共軛自同構、範數、平衡基分解、序列化 |
| lazer `polyvec.c` | `algebra::ring::polyvec::PolyVec<D>` | ✅ | 內積、範數、逐係數分解 |
| lazer `polymat.c` | `algebra::ring::polymat::PolyMat<D>` | ✅ | 行主序,矩陣-向量乘 |
| lazer `ntt.c` + `moduli.h` + hexl | `algebra::ntt` | ✅ | 見下「NTT 設計」 |
| lazer `rng.c` | `algebra::rng` | ✅ | SHAKE128(默認,比特兼容)與 AES256-CTR 兩後端,(seed, dom) 域分離 |
| lazer `shake128.c` | `sha3` crate | ✅ | 以 NIST 空消息 KAT 驗證 |
| lazer `aes256ctr*.c` | `aes` + `ctr` crate | ✅ | IV = dom_le‖0,`Ctr64LE`;與 C 的逐字節對拍待做(見「已知偏差」) |
| lazer `urandom.c` | `algebra::sampling::uniform` | ✅ | 位截斷拒絕採樣 |
| lazer `brandom.c` | `algebra::sampling::binomial` | ✅ | CBD(k),比特佈局與 C 一致 |
| lazer `grandom.c` | `algebra::sampling::gaussian` | ✅ | CDF155 常數逐字移植;FALCON 式雙峰 + BerExp |
| lazer `rejection.c` | `algebra::sampling::rejection` | ✅ | rej_standard / rej_bimodal,對數域 f64(見「已知偏差」) |
| lazer `int.c`/`intvec.c`/`intmat.c` | 暫未移植 | ⬜ | 大整數層;地基層用 i64/i128 已覆蓋 q<2^32 的需求,lazer 大模數環需要時再引入(或用 crypto-bigint) |
| lazer `coder.c`/`dcompress.c` | 暫未移植 | ⬜ | 熵編碼與 Δ 壓縮,證明序列化層 |

## 證明系統層映射(第二輪完成)

參照實現為 [lattice-dogs/labrador](https://github.com/lattice-dogs/labrador)
(LaBRADOR 核心 + Chihuahua/Dachshund/Greyhound 前端,README 稱其為
"all the dogs")。

| C 來源 | Rust 位置 | 狀態 | 備註 |
|---|---|---|---|
| `data.h` 常數(TAU1/TAU2/T/SLACK/LOGDELTA/LIFTS) | `labrador::params` | ✅ | 32 個 ±1、8 個 ±2、算子範數界 14、LIFTS=4 |
| `poly.c` `challenge()` + `poly_opnorm` | `algebra::sampling::challenge` | ✅ | inside-out Fisher–Yates + 複嵌入算子範數拒絕(Python 模擬接受率 ~12%) |
| `poly.c` `quarternary` | 同上 | ✅ | 聚合挑戰:係數均勻 {−2,−1,0,1} |
| `labrador.c` `sis_secure` | `labrador::params::sis_secure` | ✅ | MSIS 根 Hermite 因子估計,公式逐字 |
| `labrador.c` `init_proof` 參數搜索 | `labrador::params::select_params` | ✅ | 秩重組(nu 標記/分組/切塊)、f/b、fu/bu、fg/bg、κ/κ1 全部公式逐字移植;quadratic∈{0,1},C 的 quadratic==2(σ₋₁ 見證加倍,Dachshund 用)超出本輪範圍 |
| `labrador.c` `commit`/`qugarbage` | `labrador::round::commit` | ✅ | 內 Ajtai t_i、分解 t̂、外承諾 u1 = B t̂ + C ĝ |
| `jlproj.c` + `project` | `labrador::jl` | ✅ | ±1 位矩陣流式生成、拒絕循環、256 標量坍縮 + σ₋₁ 共軛;C 用 AES128 生成矩陣,Rust 用 SHAKE(見偏差) |
| `labrador.c` `lift_aggregate_zqcnst` | `labrador::round`(lifts) | ✅ | 4 次平行 Z_q 約束提升;bb 常數項置零發送,quarternary α 摺疊 |
| `labrador.c` `aggregate` | `labrador::round::aggregate_core` | ✅ | α/β/γ/δ 摺疊 u1/u2 開啟、Az=Σcᵢtᵢ、垃圾項、z⊗z 二次矩陣,偏移與縮放逐字 |
| `labrador.c` `amortize`(+tail) | `labrador::round`(amortize) | ✅ | 線性垃圾 h、u2 = D ĥ、挑戰、z = Σcᵢsᵢ;tail 模式逐 i 順序挑戰、h 明文 |
| `labrador.c` `verify` | `labrador::statement::verify` | ✅ | 全部六類檢查(範數、雙外承諾開啟、攤銷開啟、b 恆等式、⟨z,z⟩、⟨φ,z⟩)含 tail 分支 |
| `chihuahua.c` sparsecnst/prncplstmnt | `labrador::principal` | ✅ | 向量值約束(extension degree)展平為逐坐標標量約束,各自持挑戰 |
| `pack.c` composite | `labrador::composite` | ✅ | 貪心尺寸規則(proof+新witness ≥ 舊witness 即停),再試一輪 tail,末端見證明文 |
| `greyhound.c` `init_polcomctx`/`polcom_commit` | `greyhound::select_params`/`commit` | ✅ | m×n 分塊、兩層承諾、參數搜索公式逐字 |
| `greyhound.c` `polcom_eval` | `greyhound::eval` | ✅ | 部分求值 w、u2 承諾、挑戰摺疊 z、四段見證 (ẑ0,ẑ1,t̂,ŵ) |
| `greyhound.c` `polcom_reduce` | `greyhound::reduce` | ✅ | 2κ1+κ+2 條稀疏約束(u1/u2 開啟、⟨a_l,z⟩=Σcᵢwᵢ、Az=Σcᵢtᵢ、ct⟨a_r,ŵ⟩=y),xvec 求值恆等式已解析驗證 |
| `codec*.c`(熵編碼證明序列化) | 暫未移植 | ⬜ | 證明壓縮編碼層 |
| `dachshund.c`(二進制 R1CS 前端) | 暫未移植 | ⬜ | 依賴 quadratic==2 見證加倍 |

架構要點:`prove_round`/`reduce_round` 是單步歸約(承諾 → JL 投影 →
Z_q 約束提升 → 環約束聚合 → 攤銷),`Input` 枚舉統一 principal(首輪)
與 core(後續輪)兩種輸入;`composite_prove` 反覆自摺疊直至消息不再縮小,
輸出 = 各輪 RoundProof + 明文末端見證;`composite_verify` 鏈式 reduce
後對末端見證直接 `verify`。所有挑戰都由 32 字節滾動哈希態派生
(帶類型域標籤),證明者與驗證者走同一 transcript。

## NTT 設計

q = 2^32−99 之下 q−1 的 2-adicity 僅為 2,X^64+1 的負循環 NTT 不存在。
C 版做法:lazer 把大模數表示為若干 ~50 比特 NTT 友好素數之積(`moduli.h`,
SageMath 生成常數表);labrador 用 K=6 個 14 比特素數的 RNS 配合 AVX-512
16 比特通道。

Rust 版對純量 u64/u128 算術取自然選擇:兩個大素數的 CRT——
Goldilocks 素數 2^64−2^32+1(2-adicity 32)加上從 2^62 向下確定性搜索的
形如 m·2^20+1 的最大素數。CRT 模 P > 2^113,覆蓋 d ≤ 4096 時的最大乘積
係數 d·(q/2)² < 2^74。與 C 版哲學性差異:不信任硬編碼魔數——素數在首次
使用時經確定性 Miller–Rabin 驗證,單位根在運行時派生並驗序,杜絕 C 世界
「常數表悄悄壞掉」一類問題。變換本體是標準 Longa–Naehrig 位反轉孿生因子
負循環 NTT(正向 CT、逆向 GS),表按 (度數, 素數) 全局緩存。

正確性驗證:對 d ∈ {2, 64, 128, 512} 與滿幅 ±2^31 輸入,與 i128 schoolbook
負循環乘法對拍;另有環公理測試(a·X^64 = −a、共軛自同構的乘法性、
⟨a,b⟩ = ct(σ(a)·b) 等 LaBRADOR 所依賴的恆等式)。

舊 `PolynomialRing` 的乘法在雙方都是滿度(64 係數)時透明改走 NTT 引擎,
短操作數保持原 schoolbook 語義(結果長度 la+lb−1),對 labrador crate
零 API 變化。

## 已知偏差(有意為之,待辦已記錄)

1. **rejection.c 的 MPFR 128 位精度**:C 用 MPFR 做 128 位超越函數再與
   128 位定點 u 比較;Rust 版在對數域用 f64(相對誤差 ~2^−52),接受概率
   的擾動統計上可忽略,但與 C 不逐比特一致。若需嚴格一致,後續可引入
   軟浮點或 `rug`(MPFR 綁定,但破壞純 Rust)。
2. **grandom.c 的函數級 static 緩存**:C 在 `gaussian155` 裡用 static 變量
   跨調用緩存符號位,狀態跨採樣調用洩漏;Rust 版把緩存放進採樣器實例、
   每次調用重置——更乾淨且線程安全,但輸出序列與 C 不同(分布相同)。
3. **AES256-CTR 計數器佈局**:按 C 的 {dom, 0} 兩個 LE u64 佈局用
   `Ctr64LE` 復刻,尚未與 C 輸出逐字節對拍(C 需要 AES-NI 硬件,本輪
   在 arm64 沙盒無法運行 C 版)。SHAKE128 後端(默認)無此問題。
4. **參數生成**:lazer 用 SageMath 腳本離線生成 moduli/params 頭文件;
   Rust 版原則是常數運行時派生+自驗證,消滅代碼生成步驟。
5. **證明系統 transcript 不與 C 逐字節兼容**:哈希態 32 字節(C:16)、
   吸收格式為標準 LE 序列化、挑戰派生帶域標籤、JL 矩陣與承諾密鑰從
   SHAKE 派生(C:AES128 + 全局 comkey 重疊切片——Rust 用互不相交的
   域分離密鑰 A/B/C/D/W,消除密鑰重用相關性)。協議結構與所有公式
   與 C 一致,但兩者的證明不可互驗。
6. **約束表示**:C 的 extension-degree 向量值約束在 Rust 中展平為
   逐坐標標量約束(每條自持 quarternary 挑戰),數學等價。
7. **composite 停止規則**:C 按熵估計的 KB 尺寸貪心;Rust 按環元素
   計數貪心——形狀相同,切換點可能差一輪。

## 下一步(按依賴序)

1. 證明序列化與壓縮編碼(`codec*.c` 的對應物):RoundProof/CompositeProof
   目前是內存結構,需要緊湊字節編碼(z 的熵編碼、u 的位打包)。
2. quadratic==2 模式(σ₋₁ 見證加倍 + Z_q→R_q 提升向量),隨後
   Dachshund 二進制 R1CS 前端。
3. 舊 `labrador::legacy` 模塊在新核心穩定後刪除;`gaussian_generator.rs`
   的 JL ad-hoc 實現由 `labrador::jl` 取代。
4. lazer LNP 工具箱層(`abdlop.c`, `lnp-quad*.c`, `lnp-tbox.c`)——
   需要大模數環,屆時引入多素數 CRT 環描述符(`polyring_t` 的對應物)。
5. 性能:comkey 批量展開(現為逐項緩存)、amortize/commit 的 rayon 並行、
   JL 投影位運算向量化、NEON/AVX2 feature flag、Montgomery 域駐留。
6. 參數複核:sis_secure 是與 C 相同的粗估計(root Hermite 因子),
   正式安全聲明前應對照 lattice-estimator 重新標定。

## 驗證基線

兩層驗證,`scripts/` 下的參照模型獨立於 Rust 代碼可隨時重跑:

- `scripts/verify_ntt.py` — Rust NTT 引擎的 1:1 Python 孿生(同素數搜索、
  同單位根派生、同蝶形結構、同 Garner CRT),對 d ∈ {2,64,128,512} 與
  滿幅 ±2^31 輸入和整數 schoolbook 精確對拍,並驗證 CRT 界。
- `scripts/verify_gaussian.py [lazer路徑]` — CDF155 表與 C 源逐對相等
  (21/21),且以 60 位精度從 σ=1.55 半高斯第一性重推導,偏差 0 ulp;
  exp_small 核在設計區間內相對誤差 < 2.2e-16。
- `scripts/verify_samplers.py` — 採樣器的字節級可執行規範:在
  hashlib.shake_128(與 sha3 crate 同為標準 Keccak,前綴穩定)上復現
  Rust 代碼的精確字節消費順序,統計矩全部落在單元測試容差內
  (var 152.96/153.76、CBD var 0.9916、拒絕採樣接受率 0.3641 vs 1/e)。
  由於種子與流一致,這些就是 `cargo test` 將復現的數值。
- `cargo test -p algebra`(在有工具鏈的機器/CI 上運行):40+ 單元測試,
  含 NTT×schoolbook 對拍、環公理(a·X^64=−a、σ 自同構乘法性、
  ⟨a,b⟩=ct(σ(a)b))、分解重組往返、序列化往返、NIST SHAKE128 KAT、
  Miller–Rabin 已知值、域分離與確定性、採樣器矩與接受率檢驗。
- `cargo test -p labrador`:挑戰多重集與算子範數、JL 投影常數項恆等式
  ct(⟨σ₋₁(Σαπ), s⟩)=⟨α,p⟩、密鑰坍縮轉置恆等式 ⟨collapse(α,K),s⟩=⟨α,Ks⟩、
  單輪 prove/reduce 收斂(prover 與 verifier 的 h/φ/b/c 逐項相等)、
  tail 輪、兩輪鏈、quadratic 約束兩輪鏈、composite 端到端 + 三種篡改拒絕
  (見證係數、提升多項式 bb、聲明 b)。
- `cargo test -p greyhound`:求值恆等式、commit→eval→reduce→
  principal_verify 未壓縮往返、錯誤求值拒絕、壓縮鏈端到端 + 篡改 y 拒絕。
- 兩輪均在無 rustc 的沙盒(arm64,網絡白名單)中開發,編譯與 clippy/fmt
  由本地 `cargo test` 與倉庫 CI(clippy -D warnings + fmt --check)把關;
  挑戰採樣器另以 Python 精確孿生驗證(多重集、範數、~12% 接受率),
  Greyhound 的 xvec/w 重組恆等式經解析推導核驗。

# 計算条件の詳細

`classical-mc-simple` で用いているモデル・アルゴリズム・観測量の定義をまとめる。

---

## 1. 入力ファイルと役割

| ファイル | 役割 | 主なキー |
|----------|------|----------|
| **param.def** | シミュレーション共通パラメータ | Burn_in, Total_Step, Sample, num_temp, Ini_T, Delta_T, lambda, H, spin_dim |
| **lattice.def** | 格子の形とサイズ | L_x, L_y, L_z, orb_num |
| **interaction.def** | 結合の定義（ボンドリスト） | 行: `i_orb  j_orb  dx  dy  dz  J` |

---

## 2. 格子

- **サイト数**: `All_N = L_x × L_y × L_z × orb_num`
- **周期境界**: すべての方向で周期境界条件（`lattice.c` の `(jx = (ix+dx)%L_x + ...)` 等）
- **オービタル**: `orb_num = 1` のときは 1 サイト 1 スピン。`orb_num > 1` のときは単位胞あたり `orb_num` 個のスピン。
- **サイト番号**: フラット化インデックス `all_i = orb + site_index * orb_num`、`site_index = ix + iy*L_x + iz*L_x*L_y`。

---

## 3. ハミルトニアン

コード内の符号規約は **E ∝ + J (S_i · S_j)**（`interaction.def` のコメントと一致）。

- **交換相互作用（異方性あり）**
  - \( \displaystyle
    E_{\rm ex} = \frac{1}{2} \sum_{\langle i,j \rangle} J_{ij}
    \bigl( s^x_i s^x_j + s^y_i s^y_j + \lambda\, s^z_i s^z_j \bigr)
  \)
  - 和はすべてのボンドを 1 回ずつ（1/2 は二重カウントの解消）。
- **ゼーマン項**
  - \( \displaystyle E_{\rm Z} = H \sum_i s^z_i \)
- **全体**
  - \( E = E_{\rm ex} + E_{\rm Z} \)

したがって:
- **強磁性**: `J < 0`（例: interaction.def で `J = -1.0`）
- **反強磁性**: `J > 0`
- **λ**: z 成分だけ異方性（XY は λ=0、イジングは実質 λ のみ使用で spin_dim=1）。

---

## 4. スピン自由度（spin_dim）

| spin_dim | モデル | スピン表現 |
|----------|--------|------------|
| 1 | Ising | \(s^x=s^y=0\), \(s^z = \pm 1\) |
| 2 | XY | 単位円上の \((s^x,s^y)\), \(s^z=0\) |
| 3 | Heisenberg | 単位球上の \((s^x,s^y,s^z)\) |

規格化: \(|{\boldsymbol s}|=1\)（イジングは \(s^z=\pm 1\)）。

---

## 5. モンテカルロ更新

- **提案**: 各サイトを 1 回ずつ訪問し、そのサイトのスピンの新しい向きを提案。
  - Ising: 反転 \(s^z \to -s^z\)
  - XY: 一様な角度で \((s^x,s^y)\) を再サンプル
  - Heisenberg: Marsaglia 法で単位球上一様
- **受諾**: Metropolis
  - \( P_{\rm acc} = \min\bigl(1,\, e^{-\Delta E/T}\bigr) \)
- **有効場**: 隣接和 \(\boldsymbol{h}_i = \sum_j J_{ij} (\ldots)\) を `env_sx, env_sy, env_sz` に保持し、\(\Delta E = (\boldsymbol{s}_{\rm new}-\boldsymbol{s}_{\rm old})\cdot\boldsymbol{h}_i\) で O(1) 評価。

---

## 6. Exchange MC（EXMC）

- **レプリカ**: 温度インデックス `int_T = 0 .. num_temp-1` ごとに 1 レプリカ（独立なスピン配置）。
- **温度**: \( T({\tt int\_T}) = {\tt Ini\_T} + {\tt Delta\_T} \times {\tt int\_T} \)
- **交換**: 隣接する温度スロット \((i,\, i+1)\) のレプリカを確率
  \[
  P_{\rm ex} = \min\Bigl(1,\; \exp\bigl( (1/T_i - 1/T_{i+1})(E_i - E_{i+1}) \bigr) \Bigr)
  \]
  で入れ替え。odd-even の順で試行（将来の MPI 境界交換と整合）。
- **1 ステップの流れ**: 全 `int_T` で 1 sweep の Metropolis → 隣接ペアへの交換スイープ → 観測。

---

## 7. サンプリング・統計

- **Burn_in**: 平衡化ステップ数。この間は観測しない。
- **Total_Step**: 測定ステップ数。各ステップで「全温度の MC + 交換」の後に観測。
- **Sample**: 独立サンプル数。シードを `base_seed + 8945 * int_samp` で変え、各サンプルで Burn_in → Total_Step を実行。
- **平均**: 観測量は「全サンプル・全測定ステップ」で平均した値を最終結果とする。

---

## 8. 観測量の定義

- **E_per_site**: \(\langle E \rangle / N\)（サイトあたりの平均エネルギー）
- **C_per_site**: 比熱（サイトあたり）
  \[
  C = \frac{ \langle E^2 \rangle - \langle E \rangle^2 }{ N\, T^2 }
  \]
- **M2**: 磁化の 2 乗（サイト数規格化）
  \[
  M^2 = \frac{ ( \sum_i s^x_i )^2 + ( \sum_i s^y_i )^2 + ( \sum_i s^z_i )^2 }{ N^2 }
  \]
- **acceptance**: Metropolis の受諾率（全サイト・全ステップでの受諾回数 / 試行回数）
- **exchange_acceptance**: EXMC の隣接交換の受諾率（全試行に対する受諾回数の割合）
- **overlap（オーバーラップ）**: **同じ温度**での、測定開始時（MC step = 0）のスピン配置と、step = t のスピン配置の内積のサイト平均。すなわち
  \[
  Q(T) = \frac{1}{N} \sum_i \mathbf{s}_i(T,\, \mathrm{step}=0) \cdot \mathbf{s}_i(T,\, \mathrm{step}=t)
  \]
  を各ステップで計算し、測定ステップ・サンプルで平均した値。EXMC ありでは温度スロットごとに交換で別温度の配置が入るため、step 0 の配置との一致度は低くなりがち（デコレレーションが促進される）。

---

## 9. サンプル設定例（square_L16_Ising）

| 項目 | 値 | 説明 |
|------|-----|------|
| **param.def** | | |
| Burn_in | 500 | 平衡化ステップ |
| Total_Step | 3000 | 測定ステップ |
| Sample | 2 | 独立サンプル数 |
| num_temp | 6 | レプリカ（温度）数 |
| Ini_T | 1.0 | 最低温度 |
| Delta_T | 0.3 | 温度刻み → T = 1.0, 1.3, 1.6, 1.9, 2.2, 2.5 |
| lambda | 1.0 | 異方性（イジングでは z のみ） |
| H | 0.0 | 磁場 |
| spin_dim | 1 | Ising |
| **lattice.def** | | |
| L_x, L_y | 16 | 2D 正方格子の幅 |
| L_z | 1 | 2D のため 1 |
| orb_num | 1 | 1 サイト 1 スピン → N = 256 |
| **interaction.def** | | |
| ボンド | (1,0,0), (0,1,0) | 最近接のみ |
| J | -1.0 | 強磁性 |

以上が、現在のコードで一貫して使われている計算条件の詳細である。

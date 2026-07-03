# RFpeptide workflow frame

## Release packaging

The public release package is built from a safe allowlist. It includes source code, `configs/DUMMY_TARGET.json`, `cluster_profile.example.json`, and the small files under `examples/`. It excludes local credentials, real PDB/CIF files, generated workflows, retrieved results, logs, and `cluster_profile.local.json`.

Build a local release zip:

```powershell
python .\scripts\build_release.py
```

The generated archives are written to:

```text
dist/RFpeptide-workflow-frame-v<version>.zip
dist/RFpeptide-workflow-frame-latest.zip
```

Before publishing, check that no real target configuration, password file, scratch output, or retrieved result has been manually added to git.

## Cluster Settings page

Cluster connection editing has been moved out of `Cluster Dashboard` into the dedicated `Cluster Settings` tab. Use this page to edit:

- SSH host, user, port, private key, and authentication mode.
- Password login and local encrypted password storage.
- Remote upload parent and workflow directory fallback.
- PyMOL executable path and `bjobs` user.
- Result scan depth, log scan depth, result file patterns, and extra SSH args.

`Cluster Dashboard` is now reserved for runtime operations: test connection, upload workflow, submit stage, query queue status, refresh jobs, inspect pending reasons, and kill jobs.

用于把 RFdiffusion → ProteinMPNN/relax → AfCycDesign/RMSD → PyRosetta 的 peptide binder 设计流程标准化，并通过本地 GUI 管理参数、生成脚本、提交集群任务、查看任务状态和下载结果。

## 计算与设计原理

这个 workflow 面向“短肽 / 环肽 binder 设计”。核心思路是把问题拆成四个阶段：先生成可能贴合靶点表面的骨架，再为骨架设计序列，然后用结构预测和物理打分筛掉不可靠设计。

### 总体设计逻辑

输入通常是一条靶蛋白结构链，例如 `target_chain=A`，以及一个待设计的 binder 链，例如 `binder_chain=B`。`contigs` 描述 RFdiffusion 需要保留的靶点区域和新生成的 binder 长度，例如：

```text
A9-175/0 8-16
```

这里的含义是：保留靶点 `A` 链 9–175 位残基，在其旁边生成一段长度约 8–16 aa 的新 binder。当前框架默认启用 cyclic peptide 设计，因此生成和后续评估都围绕“短肽首尾可环化、同时能稳定结合靶点”的目标展开。

`RFDiffusion` 页中的 `Hotspot residues` 用于指定希望 binder 接触的靶点残基。可填写单个热点 `A326`，或多个热点 `A67,A141`；也接受备注式写法 `A326(F)`，生成脚本时会规范化为 `A326`。留空时不会写入热点参数。生成后的命令示例：

```text
'ppi.hotspot_res=[A67,A141]'
```

### 阶段 1：RFDiffusion 生成 binder 骨架

`RFDiffusion` 负责从靶点结构出发，采样大量可能的 binder backbone。这个阶段主要探索几何形状，而不是最终氨基酸序列。

它回答的问题是：

- 目标蛋白表面附近能否放置一段合理长度的短肽骨架。
- 这段短肽是否能以合适的空间朝向接触靶点。
- 在 cyclic peptide 设置下，N/C 端是否有机会靠近并形成可环化构象。

输出是一批 backbone PDB。因为扩散采样随机性较强，所以这里采用 shard 并行：例如默认 `10 shards × 20 designs`，得到足够多的候选结构供后续筛选。

### 阶段 2：ProteinMPNN + relax 设计序列

`ProteinMPNN` 给 RFdiffusion 产生的 backbone 设计氨基酸序列。它固定或近似固定 backbone 几何，寻找更可能折叠并稳定结合的序列。

随后 relax 步骤对结构进行局部能量优化，减少明显的原子冲突和不合理构象。这个阶段的重点是：

- 为同一个 backbone 生成可表达的氨基酸序列。
- 降低局部 clash 和几何异常。
- 输出后续 AfCycDesign 和 PyRosetta 可以直接评估的复合物结构。

### 阶段 3：AfCycDesign 结构预测与 RMSD 检查

`AfCycDesign` 用 AlphaFold/ColabDesign 相关模型重新预测 target-binder 复合物，用于检查设计是否在独立预测中仍能保持原始结合模式。

这个阶段不是重新设计序列，而是验证：

- binder 是否仍预测为稳定结构。
- binder 与 target 的相对位置是否接近设计输入。
- interface PAE、pLDDT 等预测置信度是否支持该结合模式。
- 预测结构和设计结构之间的 RMSD 是否过大。

如果一个候选只在原始 backbone 中看起来合理，但预测后 binder 漂移、脱离界面或置信度低，就应在这一阶段被淘汰。

### 阶段 4：PyRosetta 物理打分和界面评估

`PyRosetta` 从 Rosetta 能量函数角度评估复合物和界面。它补充的是更偏物理/几何的判断，例如：

- 复合物总能量和界面能量。
- interface ΔG、界面面积、氢键和未满足极性原子。
- packstat、SAP、接触面积等界面质量指标。
- binder 序列的疏水性、电荷、芳香性、环化可行性等辅助指标。

这个阶段用于把 AfCycDesign 中“预测可靠”的候选进一步排序，优先保留同时满足结构预测置信度和物理打分的设计。

### 筛选原则

最终结果通常不应只看单一分数。更稳健的筛选思路是组合判断：

- `AfCycDesign`：高 binder pLDDT、低 interface PAE、合理 RMSD。
- `PyRosetta`：较好的 interface energy、packstat、ΔSASA、较少 unsat。
- `序列属性`：长度适合环化，N/C 端不过度 bulky，整体电荷和疏水性不过分极端。
- `人工检查`：用 PyMOL/ChimeraX 检查 top candidates 的界面、环化距离和明显 clash。

因此，本项目的 GUI 和结果检视功能不是替代科学判断，而是把重复的生成、提交、追踪、下载和初筛流程标准化，减少手工复制脚本、改路径和漏掉任务状态的风险。

## 功能

- 编辑靶点参数：target、pilot、PDB、chain、contig、scratch 路径等。
- 在 `Project` 页填写基本信息后，一键同步各阶段的项目派生路径和 Results/Scratch 上下文。
- 编辑各阶段参数：队列、环境名、脚本路径、shard 数、recycle/model 数等。
- 一键生成四阶段 LSF 提交脚本。
- 通过 SSH/SCP 上传 workflow 到集群。
- 通过 `bsub` 提交指定阶段。
- 通过 `bjobs` 实时刷新任务状态。
- 通过 `bkill` 删除单个或批量任务。
- 查看每次 SSH/LSF 命令的完整 stdout/stderr 原始输出。
- 将 GUI 操作和远程返回写入本地 debug 日志，方便复制报错信息。
- 扫描 scratch 工作目录，按 target / pilot / stage / shard 检视结果文件。
- 在结果列表中搜索、排序，并下载选中或筛选后的文件。
- 独立的 Scratch Safety 页面：把旧日期 scratch 用户目录整体移动到最新日期目录，减少被清理风险。

## PDB 预处理

GUI 的第一个页签是 `PDB Preprocess`，用于在正式生成 workflow 前处理靶点输入结构：

- 顶部 `Local PDB` 可直接选择本地原始 PDB，并用 `Preview/Clean` 打开预处理页，避免反复手写路径。
- `Experimental RCSB PDB fetch` 可输入四位 PDB ID，例如 `8BAV`，从 RCSB 下载 legacy PDB 格式文件到 `inputs_preprocessed/rcsb/`。如果某个超大结构只提供 mmCIF，需要手动下载并转换为 PDB 后再使用当前清理功能。
- `Preview structure` 会读取蛋白质残基，按 chain 和 residue number 显示序列，并在右侧表格列出 `chain / resseq / insertion code / resname / one-letter aa / atom count`。
- `Structure components` 会以易读分区显示 PDB 标题、各链的 molecule 描述、PDB/UniProt 数据库引用、蛋白质/非蛋白质残基数量、非标准成分及其 `HETNAM` 描述、生物学 assembly 和 alternate-location 原子统计。
- `Chains to keep` 支持 Ctrl/Shift 多选；`Select all`、`Protein chains` 和 `Clear` 用于快速调整。`Apply cleanup / keep chains` 只保留选中的链。
- 清理后的 PDB 会重建保留链的 `COMPND` 描述并保留可用的 `DBREF` 信息，因此再次预览时不会丢失 chain description。
- `Rename retained chains from A` 会按原文件中的链顺序，把保留链依次重命名为 `A`、`B`、`C`……；完成后的状态栏会显示旧链到新链的映射。启用后需要确认 `target_chain`、RFdiffusion `contigs` 和 `Hotspot residues` 都使用新的链编号。
- `Remove non-protein components` 会移除水、离子、小分子配体等非蛋白质 `HETATM`/残基；`MSE` 会按蛋白质残基识别为 `M`。
- `Renumber atom serials` 只重排 atom serial，通常安全，不改变 RFdiffusion contig 使用的残基编号。
- `Renumber residues continuously by chain` 会把每条链的残基重新编号为连续编号；这会改变 RFdiffusion 的 `contigs` 所引用的 residue number，因此只有在你明确要清理不连续编号时再启用。
- 选择或 fetch 新结构时会在 `inputs_preprocessed/originals/` 建立原始快照；`Restore original PDB` 可一键恢复当前预览和 workflow 输入。
- `Apply cleanup / keep chains` 写出清理后的 PDB；勾选 `Set cleaned PDB as workflow input` 后，生成 workflow 时会自动把这个清理后的 PDB 打包到 `inputs/` 中。
- 所有主要页面现在使用统一的横向和纵向滚动容器；表格内部仍保留各自的滚动条。

## 文件结构

- `configs/PGLYRP1.json`：靶点 workflow 配置示例。
- `generate_workflow.py`：根据 JSON 配置生成四阶段提交脚本。
- `workflow_gui.py`：桌面 GUI，包括参数编辑和集群仪表盘。
- `pdb_preprocess.py`：本地 PDB 预览、序列编号查看和清理工具。
- `cluster_ops.py`：SSH/SCP/LSF 操作封装。
- `cluster_profile.example.json`：集群连接配置示例。
- `environment.yml`：推荐的 Conda/Micromamba 环境配置。
- `requirements.txt`：pip 依赖清单。
- `launch_gui.ps1`：Windows PowerShell 启动脚本。
- `generated_examples/`：示例输出目录。
- `credentials.local.json`：使用 `Remember securely` 后生成的本地加密密码文件，默认不提交。
- `logs/workflow_gui_debug.log`：运行时生成的 debug 日志，默认不提交。

## 启动 GUI

在仓库根目录运行：

```powershell
.\launch_gui.ps1
```

如果你已有 Python 环境，也可以直接运行：

```bash
python workflow_gui.py
```

## Python 环境

推荐用 Conda 或 Micromamba 创建本地 GUI 环境：

```powershell
conda env create -f environment.yml
conda activate rfpeptide-workflow-frame
```

如果使用 Micromamba：

```powershell
micromamba create -f environment.yml
micromamba activate rfpeptide-workflow-frame
```

然后启动：

```powershell
python workflow_gui.py
```

### 密码登录依赖

如果使用 Xshell 同样的账号密码登录方式，GUI 需要 Python 包 `paramiko`。安装到运行 GUI 的 Python 环境：

```powershell
python -m pip install -r requirements.txt
```

如果你用的是 Codex bundled Python，可以用完整路径安装，例如：

```powershell
C:\Users\usayz\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe -m pip install -r requirements.txt
```

默认情况下，密码只在 GUI 运行时保存在内存里，不会写入 `cluster_profile.local.json`。

如果勾选 `Remember securely`，密码会使用 Windows DPAPI 加密后保存到：

```text
credentials.local.json
```

这个文件已加入 `.gitignore`。DPAPI 密文通常只能由当前 Windows 用户在当前机器上解密；不要把它当成可跨机器迁移的密码文件。

## 配置集群连接

复制示例配置：

```powershell
Copy-Item .\cluster_profile.example.json .\cluster_profile.local.json
```

编辑 `cluster_profile.local.json`：

```json
{
  "host": "cluster.example.edu",
  "user": "bme-yaozm",
  "port": 22,
  "auth_method": "password",
  "remember_password": false,
  "ssh_key": "C:/Users/usayz/.ssh/id_ed25519",
  "remote_workflow_dir": "$HOME/RFpeptide-workflow-frame/PGLYRP1_workflow",
  "remote_upload_parent": "$HOME/RFpeptide-workflow-frame",
  "scratch_scan_roots": [
    "/scratch/2026-05-24/bme-yaozm/PGLYRP1_test"
  ],
  "scratch_scan_max_depth": 8,
  "local_download_dir": "downloads/PGLYRP1",
  "bjobs_user": "$USER",
  "scratch_migration": {
    "scratch_root": "/scratch",
    "source_date": "2026-05-24",
    "target_date": "auto",
    "user_dir": "bme-yaozm"
  }
}
```

`cluster_profile.local.json` 已加入 `.gitignore`，不要提交账号、私钥路径或服务器信息。

## 推荐使用流程

1. 在 GUI 中打开一个配置模板，或另存为新靶点配置。
2. 在 `Project` 页填写 `target`、`pilot`、cluster user、scratch date、scratch project dir 和 cluster HOME project dir。
3. 点击 `Apply project info to all pages`。它会先根据 Cluster Settings 的 `remote_upload_parent` 推导 `<remote_upload_parent>/<target>_workflow`，再把远程输入 PDB 设置为该目录下的 `inputs/<pdb name>`，把 AfCycDesign/PyRosetta 脚本设置为 `scripts/<bundled script name>`。它还会同步 Results Browser 扫描日期/用户、Scratch Safety 来源日期/用户以及 Dashboard 远程 workflow 目录。为空的 `pilot`、scratch project dir 和 HOME project dir 会分别补为 `pilot0`、`<target>_test` 和 `$HOME/<target>`。
4. 检查 RFDiffusion 的 chain、contig、hotspot residues，以及各阶段 queue、shard、recycle/model 等计算参数；同步按钮不会修改这些计算参数。
5. 点击顶部 `Save` 保存配置，再点击 `Generate workflow` 生成本地提交脚本。
6. 在 `Cluster Dashboard` 中加载 `cluster_profile.local.json`。
7. 在 `Authentication` 中选择 `password`，输入你的集群密码；如需保存，勾选 `Remember securely`。
8. 点击 `Test connection` 确认 SSH、`bjobs`、`bsub`、`bkill`、`queueinfo` 可用。
9. 点击 `Upload workflow` 上传生成目录到集群。
10. 选择阶段并点击 `Submit stage` 提交任务。
11. 在 `Queue status` 中点击 `queueinfo` 查看全部队列，点击 `queueinfo -gpu` 查看 GPU 队列；填写队列名后点击 `queueinfo -l` 查看具体队列。
12. 勾选 `Auto refresh every 10s` 实时刷新 `bjobs` 状态。
13. 需要停止任务时，选择任务后点 `Kill selected`，或点 `Kill all visible` 批量 `bkill`。
14. 任务处于 `PEND` 时，可选中任务并点击 `View bjobs -l` 查看 LSF 返回的 pending reason。

`Upload workflow` 和 `Submit stage` 会按当前打开的 workflow 配置动态派生远程目录和提交脚本名。例如当前 `target=BPIFA1` 时，即使 `cluster_profile.local.json` 里仍残留 `PGLYRP1_workflow`，GUI 也会提交 `$HOME/RFpeptide-workflow-frame/BPIFA1_workflow` 下的 `submit_BPIFA1_...sh`。提交前日志会打印实际使用的远程 workflow 目录和阶段脚本。

生成真实项目时，`local_afcyc_script`、`local_rmsd_script`、`local_merge_script`、PyRosetta `local_script` 和 `local_merge_script` 必须指向真实计算脚本。`examples/scripts/` 下的是 release 演示用 dummy entrypoint；生成器现在会拒绝把这些 dummy 脚本打包到非 `DUMMY` 靶点。AfCyc shard 运行结束后也会检查 `results.csv` 和 `pred_pdbs/*.pdb`，缺少结果时任务会明确以错误状态退出，而不是显示 `Done` 但输出为空。
15. 查看 `Raw SSH/LSF output` 面板，确认提交、删除、队列查询、下载等命令的完整返回。
16. 在 `Results Browser` 中点击 `Scan scratch` 扫描 scratch 输出。
17. 用搜索框、下拉筛选和列标题排序定位目标文件。
18. 点击 `Scan logs/errors` 专门扫描 `.out/.err` 文件，选中后可用 `Preview log text` 预览。
19. 点击 `Download selected` 或 `Download all filtered` 下载文件。

## Scratch Safety

`Scratch Safety` 是独立隔离功能，不参与任务提交、任务删除或结果下载。

用途：把旧 scratch 日期目录下的用户文件夹整体移动到目标日期目录，避免旧日期目录被集群清理。例如：

```text
/scratch/2026-05-24/bme-yaozm
```

移动到：

```text
/scratch/2026-06-05/bme-yaozm
```

配置项：

- `scratch_root`：通常是 `/scratch`。
- `source_date`：当前旧日期目录，例如 `2026-05-24`。
- `target_date`：目标日期目录；填 `auto` 时使用集群端 `date +%F`。
- `user_dir`：你的 scratch 用户目录，例如 `bme-yaozm`。

安全规则：

- `Use current target settings` 会从当前 workflow 配置填入 `source_date` 和 `user_dir`。
- `Preview move` 只检查源目录、目标目录、目录大小和顶层文件列表，不移动数据。
- `Move now` 执行前会弹窗确认。
- 如果源目录不存在，会拒绝执行。
- 如果目标用户目录已经存在，会拒绝执行，不会合并或覆盖。
- 执行结果会完整写入 `Raw SSH/LSF output` 和 `logs/workflow_gui_debug.log`。

## 结果检视

`Results Browser` 的扫描目录由界面上的 `Scratch root`、`Date` 和 `User folder` 决定，实际扫描：

```text
<scratch root>/<date>/<user folder>
```

例如：

```text
/scratch/2026-06-09/bme-yaozm
```

点击 `Use current target` 会从当前 workflow 配置的 `scratch_date` 和 `cluster_user` 自动填入。旧版 `cluster_profile.local.json` 中的 `scratch_scan_roots` 仍作为初始值读取，但扫描时会被界面上的日期/用户设置覆盖。

默认识别：

- `*.csv`
- `*.pdb`
- `*.out`
- `*.err`
- `*.json`
- `*.txt`

扫描结果会尽量从路径中推断：

- `target`：例如 `PGLYRP1_test` 或 `PGLYRP1_workflow`。
- `pilot`：例如 `pilot0`。
- `stage`：`RFDiffusion`、`ProteinMPNN`、`AfCycDesign`、`PyRosetta`、`Logs` 或 `Other`。
- `shard`：例如 `shard00`、`runlist_03` 这类路径中的编号。

扫描时 `Status` 会显示当前正在扫描的根目录；扫描结束后会显示文件数量。如果结果为 0，日志窗口也会记录未扫描到文件的根目录，便于判断是路径不对还是确实没有输出。

RFDiffusion 的 PDB 会进一步细分：

- 正式设计：`RFDiffusion`，例如 `diffused_binder_cyclic_mcl1_shard00_0.pdb`。
- `pX0` 扩散轨迹：`RFDiffusion pX0 traj`，例如 `diffused_binder_cyclic_mcl1_shard00_0_pX0_traj.pdb`。
- `Xt-1` 扩散轨迹：`RFDiffusion Xt-1 traj`，例如 `diffused_binder_cyclic_mcl1_shard00_0_Xt-1_traj.pdb`。

下载时会自动整理到：

```text
downloads/<target>/<pilot>/<stage>/<shard>/
```

在 GUI 的 `Results Browser` 中，`Download selected` 和 `Download all filtered` 会优先下载到当前 workflow 目录下的新建独立目录：

```text
generated_examples/<target>_workflow/retrieved_results/<stage>/<target>/<pilot>/<shard>/
```

因此 PDB 和 CSV 都会先按 `stage` 分组，再保留 `target/pilot/shard`，避免不同阶段或 shard 的同名文件互相覆盖。CSV 打分结果通常来自：

- `AfCycDesign`：`afcyc_out/<pilot>/merged/*.csv`
- `PyRosetta`：`pyrosetta_scores/<pilot>/merged/*.csv`

选中一个 `.csv` 后点击 `Preview CSV` 会自动下载并在内置 CSV 表格窗口中打开。CSV 预览采用分页读取，默认每页 200 行，不会一次性把大型打分表全部载入内存；左侧 `Columns` 面板可搜索列名，点击列名会把表格横向滚动到对应列。选中一个 `.pdb` 后点击 `Open PDB in PyMOL` 会自动下载并调用本机 PyMOL 打开。

## 初步数据处理

`Data Processing` 页只扫描 merged CSV，例如 `results_merged_<target>_<pilot>.csv` 和 `pyrosetta_scores_merged_<target>_<pilot>.csv`。选中 CSV 后可以：

- `Preview selected CSV`：用分页 CSV 浏览器查看大表。
- `Load for analysis`：检测可数值化的指标列。
- `Plot histogram + scatter`：对选中指标同时画直方图和按 CSV 顺序排列的散点图，并默认标记当前显示范围内的 percentile 分位线。
- `AfCyc defaults` / `PyRosetta defaults`：自动选择常用打分项和方向，例如 AfCycDesign 的低 iPAE/RMSD、PyRosetta 的低 interface_dG 或高 CMS。
- `Direction`：设置候选条目是按 `higher` 还是 `lower` 方向筛选；PyRosetta 的 `interface_dG`、`SAP`、unsat 等通常用 `lower`。
- `Hist min` / `Hist max`：手动指定直方图和散点图的显示上下限；百分位数和默认候选阈值会基于这个范围内的数据重新计算。
- `Auto robust range 1-99% when min/max blank`：上下限留空时只用 1–99% 范围绘图，减少极端异常值对尺度的扭曲；候选筛选仍基于完整数据。
- `Custom threshold`：填写后按自定义阈值列举候选；为空时按 percentile 阈值列举。

候选表会按所选方向列出通过阈值的 top 条目，包括 CSV 行号、指标值、自动推断的候选 ID 和若干摘要字段。`Data Processing` 的 CSV 列表有独立的 `Search`、`Target`、`Pilot`、`Stage` 和 `Shard` 筛选状态，不会影响 `Results Browser`。

### Hit 筛选

`Data Processing` 页下方的 `Hit screening from AfCycDesign + PyRosetta merged CSV` 用于把 AfCycDesign 和 PyRosetta 的 merged CSV 合并筛选。典型用法：

1. 在 `Merged CSV scan` 中扫描 scratch，并分别选中 AfCycDesign 和 PyRosetta 的 merged CSV。
2. 在 hit screening 的 `CSV inputs` 中分别点击 `Use selected scanned CSV`，程序会自动下载远端 CSV 并填入本地路径；也可以用 `Browse` 选择本地 CSV。
3. 使用默认阈值或手动调整列名/阈值。默认参考 RFpeptides 论文中的筛选逻辑：低 AfCycDesign iPAE、低设计模型与 AfCycDesign 预测 RMSD、低 Rosetta/PyRosetta interface dG、低 SAP、高 CMS。`binder pLDDT` 阈值留空时不参与筛选。
4. 点击 `Run hit screening`。命中结果按 interface dG、iPAE、RMSD、SAP 和 CMS 排序。
5. 点击 `Save hit CSV` 导出 `generated_examples/<target>_workflow/retrieved_hits/<timestamp>_csv/hit_screening_filtered.csv`。
6. 选中 hit 后点击 `Download selected hit structures`，或点击 `Download all hit structures`，下载对应 AfCycDesign 预测结构和 MPNN 输入结构。

下载结构会保存到：

```text
generated_examples/<target>_workflow/retrieved_hits/AfCycDesign_predicted/
generated_examples/<target>_workflow/retrieved_hits/MPNN_input/
```

### 任务日志和报错文件

当前生成的 LSF `.out/.err` 不是统一写在上传后的 workflow 脚本目录中，而是按各阶段实际工作目录写到 scratch：

- `RFDiffusion`：`/scratch/<date>/<user>/<project>/` 下的 `*.%J.out` 和 `*.%J.err`。
- `ProteinMPNN`：`/scratch/<date>/<user>/<project>/mpnn_relax*_out/<pilot>/logs/`。
- `AfCycDesign`：`/scratch/<date>/<user>/<project>/afcyc_out/<pilot>/logs/`。
- `PyRosetta`：`/scratch/<date>/<user>/<project>/pyrosetta_scores/<pilot>/logs/`。

`Results Browser` 的 `Scan logs/errors` 会从 `cluster_profile.local.json` 的 `job_log_scan_roots` 扫描；如果未配置，则复用 `scratch_scan_roots`。扫描到的 `.out/.err` 文件可以用 `Preview log text` 自动下载并预览。

### PyRosetta shard 数量检查

PyRosetta 阶段会在提交时扫描：

```text
/scratch/<date>/<user>/<project>/mpnn_relax*_out/<pilot>/*.pdb
```

然后按当前这一刻扫描到的 PDB 文件生成 runlist。实际提交的 shard 数量等于“非空 runlist 数量”，因此如果上游 ProteinMPNN/relax 输出还在增长，或者输入 PDB 数量小于配置的 `pyrosetta.n_shards`，实际提交数可能少于期望值。

为避免在上游仍在写文件时误提交，PyRosetta 提交脚本现在会：

- 先统计输入 PDB 数量。
- 等待 `INPUT_STABILITY_SECONDS`，默认 30 秒。
- 再次统计输入 PDB 数量。
- 如果数量变化，直接停止提交，并提示等待上游完成后重新提交。
- 如果输入 PDB 数量少于 `N_SHARDS`，会明确警告只能提交对应数量的非空 shard。

如果确认上游已经完成但不想等待稳定性检查，可在集群手动运行时临时设置：

```bash
WAIT_FOR_STABLE_INPUT=0 bash 4.PyRosetta/submit_pyrosetta_shards.sh
```

### GPU 资源请求

RFDiffusion 和 AfCycDesign 的 GPU shard 默认使用同一套 LSF 资源写法：

```text
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -gpu "num=1:aff=no"
```

这里使用 `span[ptile=8]` 并把 GPU affinity 显式关闭为 `num=1:aff=no`，避免 LSF 产生无法满足的 host affinity 约束。RFDiffusion 页和 AfCycDesign 页都可以分别调整 `gpu_req`、`gpu_ncpu` 和 `gpu_span`；当前默认值为 `gpu_req=num=1:aff=no`、`gpu_ncpu=8`、`gpu_span=span[ptile=8]`。生成后的 `#BSUB` 资源行会直接写成这些值，不再被登录环境里的同名变量覆盖。

如果 PyMOL 不在系统 PATH 中，可在 `cluster_profile.local.json` 中加入：

```json
{
  "pymol_executable": "C:/Program Files/PyMOL/PyMOLWin.exe"
}
```

## Debug 日志

GUI 会把操作日志和远程命令原始返回写入：

```text
logs/workflow_gui_debug.log
```

如果需要排错，把这个文件中对应时间段的内容发出来即可。它会包含：

- 操作名称。
- 实际执行的本地 SSH/SCP 命令。
- 退出码。
- stdout。
- stderr。

## 默认 shard 配置

当前默认值沿用 `Test3_PGLYRP1`：

- `RFDiffusion`：10 shards × 20 designs。
- `ProteinMPNN`：20 shards。
- `AfCycDesign`：15 shards。
- `PyRosetta`：50 shards。

## 命令行生成

不使用 GUI 时可运行：

```bash
python generate_workflow.py --force
```

指定新配置：

```bash
python generate_workflow.py --config configs/NEW_TARGET.json --force
```

## 注意事项

- GUI 不依赖 Xshell。密码模式通过 `paramiko` 连接集群；SSH key 模式通过标准 `ssh`、`scp` 调用集群。Xshell 仍可用于人工排错。
- 密码不会写入 cluster profile；如启用记住密码，只会写入 DPAPI 加密后的 `credentials.local.json`。
- `bkill` 是破坏性操作，GUI 会要求确认后再执行。
- `Scratch Safety` 的 `Move now` 也是破坏性操作，必须先 `Preview move` 并确认路径正确。
- `Download results` 仍保留为固定文件下载；更推荐使用 `Results Browser` 扫描、筛选并下载。

## 上传内容与远程结构

点击 `Upload workflow` 时，GUI 会先根据当前参数重新生成本地 workflow 包，然后把整个目标目录递归上传到集群。默认配置中，本地目录是：

```text
generated_examples/PGLYRP1_workflow
```

如果 `cluster_profile.local.json` 中配置为：

```json
{
  "remote_upload_parent": "$HOME/RFpeptide-workflow-frame",
  "remote_workflow_dir": "$HOME/RFpeptide-workflow-frame/PGLYRP1_workflow"
}
```

GUI 会优先使用 `remote_upload_parent` 作为父目录，并用当前 `target` 生成实际远程目录名，例如 `BPIFA1_workflow`。因此切换靶点后不需要手动改 `remote_workflow_dir`；如果点击 `Save` 保存 cluster profile，当前靶点派生出的 `remote_workflow_dir` 会写回配置文件。

上传后的远程目录就是：

```text
$HOME/RFpeptide-workflow-frame/PGLYRP1_workflow
```

当前生成器会把 workflow 自身需要的提交脚本、项目本地 Python 脚本和输入结构一起打包：

```text
PGLYRP1_workflow/
  inputs/
    1YCK_clean.pdb
  scripts/
    afcyc_predict_batch.py
    rmsd_from_afcyc.py
    merge_afcyc_csvs.py
    PyRosetta_fullScoring_v4_debug.py
    merge_pyrosetta_csvs.py
  1.RFDiffusion/
    submit_PGLYRP1_rfdiffusion_pilot0.sh
  2.ProteinMPNN/
    submit_PGLYRP1_mpnn_relax4_shards.sh
  3.AfCycDesign/
    submit_afcyc_shards_integrated.sh
  4.PyRosetta/
    submit_pyrosetta_shards.sh
```

`RFDiffusion`、`AfCycDesign` 和 `PyRosetta` 的提交脚本会从自身位置计算：

```bash
WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
```

因此这些脚本不再依赖旧的固定路径 `$HOME/Test3_PGLYRP1/...` 来寻找打包进 workflow 的 PDB、AfCycDesign 脚本、RMSD 脚本、PyRosetta 脚本和 merge 脚本。

不在上传包内、仍要求集群环境预先存在的内容包括：

- `RFDiffusion` 主程序目录，例如 `$HOME/RFdiffusion`。
- `ProteinMPNN/relax` 主程序，例如 `$HOME/dl_binder_design/mpnn_fr/dl_interface_design.py`。
- AlphaFold/ColabDesign/Rosetta/PyRosetta 相关模型、权重、数据库和已配置环境。
- `micromamba`、CUDA/module、LSF 的 `bsub`、`bjobs`、`bkill` 等集群基础设施。

如果以后希望 ProteinMPNN 阶段也完全脱离 `$HOME/dl_binder_design/...` 的固定路径，需要把该仓库或对应脚本及其依赖也纳入打包策略；当前版本把它视为集群已安装的外部计算环境。

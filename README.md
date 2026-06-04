# Binder-design workflow template

这个目录把 `Test3_PGLYRP1` 中已经跑通的集群提交流程抽成一个可复用框架。

## 文件

- `configs/PGLYRP1.json`：一个靶点的全部可变参数。
- `generate_workflow.py`：根据配置生成四阶段提交脚本。
- `generated_examples/`：默认输出目录，不建议手工维护其中内容。

## 默认 shard 配置

当前默认值沿用 `Test3_PGLYRP1`：

- `RFDiffusion`：10 shards × 20 designs。
- `ProteinMPNN`：20 shards。
- `AfCycDesign`：15 shards。
- `PyRosetta`：50 shards。

## 生成脚本

### GUI 方式

在 Windows PowerShell 中运行：

```powershell
.\Qiming_input\workflow_template\launch_gui.ps1
```

GUI 可以完成：

- 打开/保存 JSON 配置。
- 编辑 target、路径、chain、contig、环境、队列、shard 数等参数。
- 一键生成四阶段提交脚本。

### 命令行方式

在本地项目目录运行：

```bash
python3 Qiming_input/workflow_template/generate_workflow.py --force
```

会生成：

```text
Qiming_input/workflow_template/generated_examples/PGLYRP1_workflow/
```

把生成目录上传到集群后，按阶段提交：

```bash
bash 1.RFDiffusion/submit_PGLYRP1_rfdiffusion_pilot0.sh
bash 2.ProteinMPNN/submit_PGLYRP1_mpnn_relax4_shards.sh
bash 3.AfCycDesign/submit_afcyc_shards_integrated.sh
bash 4.PyRosetta/submit_pyrosetta_shards.sh
```

## 新靶点使用方式

复制 `configs/PGLYRP1.json`，修改这些字段：

- `target`
- `project_dir_name`
- `home_project_dir`
- `input_pdb`
- `target_chain`
- `binder_chain`
- `contigs`
- `scratch_date`
- `afcyc.*_script`
- `pyrosetta.script`
- `pyrosetta.merge_script`

然后指定配置生成：

```bash
python3 Qiming_input/workflow_template/generate_workflow.py --config Qiming_input/workflow_template/configs/NEW_TARGET.json --force
```

如果你的系统只提供 `python` 命令，也可以把上面的 `python3` 换成 `python`。

## 注意

- 生成器只负责生成提交脚本；`afcyc_predict_batch.py`、`rmsd_from_afcyc.py`、`merge_afcyc_csvs.py`、`PyRosetta_fullScoring_v4_debug.py` 等实际计算脚本仍需要放在配置指定路径。
- PyRosetta 生成脚本统一使用 `PyRosetta_fullScoring_v4_debug.py`，避免旧脚本中外层和 job 内部路径不一致的问题。

# AquaMateR

> Hatchery-aware genomic mate allocation for aquaculture breeding.

AquaMateR 是一个面向水产育种的基因组选配 R 包原型。它的目标不是简单预测 GEBV，也不是把畜牧选配算法改名后直接使用；它面向孵化场的真实决策：**谁和谁配、采用何种繁殖模式、每个组合投放多少受精卵、预计保留多少苗种，以及怎样约束实际家系贡献。**

当前版本为 **0.5.0**。它将“交配组合数”与“后代配额”分开，家系和亲本限制按预计保留苗种的遗传贡献计算，并将基于基因型分离的后代潜力、严格 OCS 与孵化场配额耦合，适合高繁殖力、家系贡献容易失衡的水产物种。

## 从 GitHub 安装

```r
install.packages("remotes")  # 只需第一次安装时运行
remotes::install_github("liyin2/AquaMateR")
library(AquaMateR)
```

安装后可运行最小示例：

```r
result <- run_aquamater_demo()
result$summary
```

> AquaMateR 当前为研究型原型软件。请在具体物种和育种群体中校准繁殖力、存活、繁殖成功偏斜与亲子鉴定策略后，再将结果用于生产决策。

## 快速运行

```r
# 在 R 中进入 AquaMateR 所在目录后：
# install.packages("devtools")
# devtools::load_all(".")

result <- run_aquamater_demo()
head(result$plan)
result$summary
```

## 核心流程

```text
SNP 矩阵 M
  ↓
calc_grm(): 构建基因组关系矩阵 G
  ↓
fit_gblup(): 估计 GEBV
  ↓
make_candidate_crosses(): 生成候选交配组合
  ↓
calc_cross_metrics(): 计算增益、近交风险、后代多样性
  ↓
predict_aqua_offspring(): 后代均值、近交与分离潜力代理指标
  ↓
allocate_hatchery_progeny(): 受精卵与苗种配额
  ↓
optimize_aqua_plan(): 按后代贡献约束优化
  ↓
输出可执行的孵化场繁殖批次计划
```

## 水产专用能力

- 支持一对一交配、因子交配与亲本池/群体产卵：`spawning_design`
- 支持每个组合的受精卵数量、早期存活率和目标保留苗种数
- 家系贡献按预计保留苗种计算，而非按交配组合计数
- 支持性别未知和雌雄同体物种，但默认不让未知性别自动进入双亲本池
- 支持同家系规避、亲本遗传贡献上限和家系遗传贡献上限
- 输出后代均值、近交代理、孟德尔分离潜力代理和后代潜力代理
- 提供多代策略比较的透明规划投影

## v0.2 基础工作流

```r
# 1. dat$M 和 dat$pheno 可以替换为真实 SNP、家系和 GEBV 数据
dat <- simulate_aqua_data(n = 120, m = 800, n_families = 20)
G <- calc_grm(dat$M, ids = dat$pheno$id)
fit <- fit_gblup(dat$pheno$y, G, ids = dat$pheno$id)

# 2. 形成孵化场亲本表；fecundity 和 survival 可来自物种/亲本记录
broodstock <- new_aqua_broodstock(
  id = dat$pheno$id,
  sex = dat$pheno$sex,
  family = dat$pheno$family,
  gebv = fit$gebv,
  fecundity = 20000,
  survival = 0.35
)

# 3. 生成并预测交配组合
crosses <- make_aqua_crosses(broodstock, spawning_design = "factorial")
outcomes <- predict_aqua_offspring(crosses, broodstock, G)

# 4. 分配保苗目标；再按家系/亲本的实际后代贡献约束选择批次
allocation <- allocate_hatchery_progeny(outcomes, target_selected = 12000)
plan <- optimize_aqua_plan(
  allocation,
  n_crosses = 30,
  max_parent_contribution = 0.10,
  max_family_contribution = 0.15
)
evaluate_aqua_plan(plan)
project_aqua_generations(plan, n_generations = 10)
```

## v0.3：把计划和实际繁殖结果闭环

```r
# 因子交配：同一雌鱼的有限卵量分配给多个雄鱼，而非把每条边都当作完整一窝。
factorial_crosses <- make_aqua_crosses(broodstock, spawning_design = "factorial")
factorial_outcomes <- predict_aqua_offspring(factorial_crosses, broodstock, G)
factorial_batches <- allocate_factorial_eggs(factorial_outcomes, n_sires_per_dam = 3)

# 亲本池/群体产卵：预演不均衡繁殖成功带来的贡献风险。
pool_sim <- simulate_parent_pool(broodstock, n_offspring = 5000, n_sim = 500)
pool_sim$summary

# 产后：将亲子鉴定结果与计划中的家系/亲本贡献进行审计。
# parentage 每行对应一个子代，至少含 sire 与 dam 两列。
audit <- audit_realised_parentage(plan, parentage, broodstock)
audit$parent
audit$family
```

`simulate_parent_pool()` 用繁殖成功率和雌性繁殖力生成群体产卵的实际贡献分布；它用于规划和敏感性分析，不能替代真实亲子鉴定。`audit_realised_parentage()` 的作用是将真实亲本贡献反馈到下一轮选配。

## v0.4：后代潜力与严格 OCS 基准

```r
# 后代潜力：输入 SNP 基因型和已估计的标记效应，模拟同一家系内的分离。
# 标记效应可以来自 BayesR、rrBLUP、GWAS 后验效应或其他遗传评估流程。
uc <- predict_usefulness_criterion(
  crosses = head(crosses, 20), M = dat$M,
  marker_effects = stats::rnorm(ncol(dat$M)),
  selection_proportion = 0.10, n_progeny = 500, seed = 1
)

# 严格 OCS：先获得“亲本应贡献多少”的凸二次规划解，
# 再由孵化场选配模块将其翻译为交配与苗种批次。
oc <- optimize_oc_contributions(fit$gebv, K = G / 2,
                                max_contribution = 0.10,
                                coancestry_penalty = 1)
```

`predict_usefulness_criterion()` 在未相位的 0/1/2 基因型下假设标记间独立分离；它适合作为 v0.4 的家系内选择基准。真实论文分析中应进一步引入物种特异的重组率、相位单倍型与验证性状的标记效应。`optimize_oc_contributions()` 是亲本贡献层的凸二次规划 OCS 基准，并不直接产生配对方案。

## v0.5：OCS - 后代潜力 - 孵化场配额耦合

```r
# 将 Usefulness Criterion 合并进候选交配组合。
outcomes$usefulness_criterion <- uc$usefulness_criterion[
  match(outcomes$cross_id, uc$cross_id)
]

# OCS 决定贡献，UC 决定同等贡献配额下优先配哪一对，
# 最终得到目标保苗量与投卵量。
oc <- optimize_oc_contributions(fit$gebv, K = G / 2,
                                max_contribution = 0.10)
coupled_plan <- allocate_oc_mating_plan(
  outcomes, oc, target_selected = 12000,
  objective = "usefulness_criterion"
)
evaluate_aqua_plan(coupled_plan)
```

这一步形成了 AquaMateR 的核心闭环：**OCS 决定亲本应贡献多少，Usefulness Criterion 决定优先形成哪类后代，孵化场分配器决定每个组合的投卵量与保苗量。**

## 跨物种使用：通用核心，物种校准

```r
# 模板只定义繁殖结构，不擅自假定任何物种的卵量或存活率。
fish_profile <- aqua_species_template("finfish")
shrimp_profile <- aqua_species_template("shrimp")
bivalve_profile <- aqua_species_template("bivalve")

# 使用真实育种群体参数创建配置，再将其用于模拟和孵化场约束。
profile <- new_aqua_species_profile(
  name = "target_population",
  spawning_design = "parent_pool",
  generation_structure = "overlapping",
  parentage_strategy = "required",
  calibration = list(
    female_egg_output = "estimate from hatchery records",
    early_survival = "estimate by batch and family",
    reproductive_skew = "estimate from parentage records"
  )
)
```

包的核心算法适用于高繁殖力水产物种；但卵量、存活、性比、倍性、群体产卵和亲子鉴定策略必须基于具体物种与育种群体校准。模板不是物种参数数据库。

`expected_inbreeding_proxy` 与 `mendelian_sampling_proxy` 是当前原型的比较指标。论文级应用应以经过标定的共祖/亲缘矩阵和物种特异的后代遗传方差模型替换它们。

## 当前版本定位

这是 AquaMateR 0.5.0 研究原型版，适合用于：

- 青年基金项目汇报
- R 包开发起点
- 模拟实验
- 方法学论文预研
- 真实水产数据试跑

如果要冲击高水平文章，建议下一步加入：

- NSGA-II 与 OCS 的混合优化
- 个体层面的多代基因组模拟
- 物种特异的 Usefulness Criterion 后代潜力模型
- G×E 多环境选配与同胞挑战试验模块
- 抗病、生长、存活率等多性状指数

## 如何引用

在 R 中安装后可运行：

```r
citation("AquaMateR")
```

推荐引用格式：

> Li L (2026). *AquaMateR: Hatchery-Aware Genomic Mate Allocation for Aquaculture Breeding*. R package version 0.5.0. Xiamen University. https://github.com/liyin2/AquaMateR

仓库根目录中的 `CITATION.cff` 可被 GitHub 和文献管理工具识别。

## 许可证与贡献

本仓库采用 **AquaMateR Research-Use No-Derivatives License**：允许研究、教学、评估和内部使用未经修改的软件；禁止修改、派生开发和发布修改版，除非事先获得作者书面许可。详见 [LICENSE](LICENSE)。

如需报告问题、提出功能建议或申请修改/派生许可，请使用 GitHub Issues 或联系 yinli@stu.xmu.edu.cn。
# AquaMateR

Genomic mating design for highly fecund aquatic species.

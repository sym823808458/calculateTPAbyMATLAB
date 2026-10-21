# Permanent-dipole reference correction / 永久偶极矩参考修复

Confirmed in source on 2026-10-21. The MATLAB overwrite loop was corrected on 2026-10-20.

## Error and correction

The first iteration of the old loop set the ground-state diagonal to zero. Later
iterations therefore subtracted zero instead of the original ground dipole. The
intended difference is Δμᵢ = μᵢᵢ − μ₀₀ for every diagonal, including the ground state.
Save the reference **before** changing the matrix:

```python
mu0 = mu[0, 0, :].copy()
for i in range(mu.shape[0]):
    mu[i, i, :] -= mu0
```

The `.copy()` is essential: without it, a NumPy slice can remain a view of the
matrix and be overwritten during subtraction. The equivalent MATLAB repair is:

```matlab
mu0 = mu(1,1,:);
for i = 1:size(mu,1)
    mu(i,i,:) = mu(i,i,:) - mu0;
end
```

Off-diagonal transition dipoles are unchanged. All permanent vectors being
subtracted must have consistent units, coordinates/axes/origin and electronic or
total definitions. This MATLAB script does not automatically normalize definitions.
**The reference-overwrite repair is independent of relaxed-density treatment.**

## Isolated impact comparison

Four existing molecular Multiwfn exports, each with 30 excited states, were
reprocessed with only the reference loop changed. Excitation energies, state
count, damping, grid and the other SOS expressions were held fixed.

| Case | Main-peak σ change |
| --- | ---: |
| 1 | +0.009236% |
| 2 | +0.855914% |
| 3 | +0.784894% |
| 4 | −0.884634% |

Main-peak positions did not change on the sampled grid. Some non-main-peak samples
changed by up to 10.72% relative to their original value; this was evaluated among
samples with original response at least 10% of the original peak. These four
small peak effects **cannot be generalized to all molecules**. There is no
universal correction factor: recompute affected inputs with the corrected code.
Original outputs and pre-fix source backups remain preserved privately. The
underlying molecular files and spectra are not redistributed here.

A common vector added to every permanent dipole left the corrected Python
spectrum invariant to numerical rounding (maximum 3.55e-15 GM in the original
audit); the unfixed implementation differed by 15.7054 GM. This checks reference
handling, not absolute cross-section accuracy or Python/MATLAB equivalence.
MATLAB runtime validation was not performed. SOS prefactors, state convergence,
damping and negative response values still require scientific evaluation.

## 中文说明

这是明确的程序错误：旧循环第一次迭代就清零基态参考，导致后续激发态实际减零。
正确实现必须在修改矩阵前保存基态独立副本；Python 的 `.copy()` 必须保留。
仅改变对角元，非对角跃迁偶极矩保持不变，并核对单位、坐标系与电子／总偶极矩定义。

四例仅修复此处、保持其余设置不变，主峰截面约变化 +0.009%、+0.856%、+0.785%、
−0.885%；部分非主峰采样点相对变化达 10.72%（原响应至少为原峰高 10% 的范围）。
影响需逐分子评估，不能外推所有分子，也不能用统一系数修正旧结果。历史输出保留，
受影响的数据应重新计算。这一修复与弛豫密度处理是独立问题，不构成整套 SOS
绝对截面的验证；MATLAB 尚未做运行验证。

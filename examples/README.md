# Reproducible example designs

These JSON files record the user-facing parameter values for representative
TEMPO designs. Every value is inside the same bounds enforced by the web
interface and Python adapter.

| File | Purpose |
|---|---|
| `baseline.json` | K4 reference used by `test_installation.py` |
| `fast_oscillator.json` | Shorter K3 oscillator setting |
| `slow_oscillator.json` | K5 setting with the certified 12 h⁻¹ tag |
| `strong_shutdown.json` | Stronger sRNA production than the reference |
| `weak_shutdown.json` | Weaker sRNA production and a longer output tail |

To reproduce an example in the interface, select the corresponding K setting
and enter the listed values. Export Design Summary from the interface to record
the resulting model version, parameters, evaluation, and system summary.

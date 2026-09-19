# Engineering the BM3R1–sRNA Shutdown Module

This dry-lab engineering study addressed a local timing problem in TEMPO: how can downstream expression be suppressed while the DNA state that enabled it remains in LR? We organized the work into two connected design–build–test–learn cycles, progressing from a functional shutdown mechanism to the conditions under which that mechanism remains effective. The current implementation models one PB/LR bit; integration with a complete target-count decoder remains a subsequent design step.

## Engineering Cycle 1 — From Persistent State to a Tunable Expression Window

### Design

**Design goal: preserve the DNA-state memory needed for counting while introducing separate control over the duration of downstream expression.**

The shutdown module needed to meet three requirements: preserve upstream counting fidelity, allow an initial expression window after PB→LR switching, and suppress new output production before the bit leaves LR.

We selected a regulatory-layer architecture that uses residual BM3R1 as a molecular timer. Following the PB→LR transition, BM3R1 production ceases as the PB state is depleted, while the existing BM3R1 pool declines over time. During this interval, BM3R1 continues to repress the shutdown promoter, allowing LR-driven target mRNA to support output production. As BM3R1 repression weakens, shutdown sRNA accumulates, pairs with the target mRNA, and reduces new translation.

This design uses an existing counter-associated timing process. Its intended shutdown action is suppression of new protein synthesis; it does not directly remove protein already produced.

![Delayed sRNA shutdown mechanism](assets/00_shutdown_mechanism.png)

*Circuit concept: LR-driven expression begins while residual BM3R1 delays sRNA activation; subsequent sRNA accumulation suppresses new translation.*

### Build

We translated this architecture into a mechanistic ODE model coupled to the existing oscillator–counter framework. The counter supplied the time-dependent LR fraction and BM3R1 level, while three additional states represented shutdown sRNA, target mRNA, and a normalized protein readout.

The model included BM3R1-dependent promoter repression, sRNA and mRNA turnover, sRNA–mRNA pairing, translation, and protein dilution. We used the RyhB/sodB system as a literature-based kinetic proxy. Crucially, the implemented pairing reaction consumes both target mRNA and sRNA: the sRNA equation includes the co-degradation term, so sRNA is not treated as an inexhaustible inhibitor.

We also represented BM3R1 binding to the added shutdown operators and its effect on RDF regulation. This allowed us to test whether adding the shutdown cassette altered counter behavior. Numerical implementation was checked against the analytical steady state of the same sRNA–mRNA equations; this was an implementation check, not experimental validation.

### Test

We compared the coupled system with a matched no-shutdown control and evaluated whether output activated, subsequently declined within LR, and retained a usable initial peak.

We defined the expression window using the full width at half maximum (FWHM) of **new translation flux**. A parameter set passed the shutdown criteria when:

- the translation-flux FWHM was finite and positive;
- mean translation flux over the final 20% of the LR interval was below 20% of that interval's peak;
- the peak remained at least 50% of the corresponding no-shutdown control.

These are project-defined engineering criteria. Counter fidelity was audited separately by checking for exactly one complete reversal per input cycle.

At the nominal 50-minute doubling time, shutdown reduced translation-flux FWHM from **10.61 h to 6.98 h**. Late-LR translation fell to **17.0% of peak**, while the initial peak retained **98.6%** of the control. With the modeled operator load included, **5/5 audited cycles** retained exactly one complete reversal.

We then scanned sRNA transcription strength and the shutdown repression threshold relative to the RDF threshold. **91/169 combinations** met the shutdown criteria, with passing FWHM values of **5.42–7.36 h**. Thus, the model produced a range of usable expression windows rather than a single working parameter point.

![Baseline simulation and no-shutdown control](assets/01_baseline_timecourse.png)

*The simulated time course tests whether new translation declines while LR remains active. Stable protein concentration is a separate readout.*

![Expression-window parameter scan](assets/02_tunability_heatmap.png)

*The two-parameter scan identifies combinations that satisfy both late suppression and peak retention.*

### Learn

The simulations showed that DNA-state persistence and continued output production need not have the same duration. Under the tested conditions, BM3R1 clearance delayed sRNA activation sufficiently to preserve an initial output peak, after which sRNA reduced new translation before the LR interval ended.

This established a useful division of temporal control in the modeled module: **the counter state enables expression, while the downstream shutdown layer helps determine how long substantial new translation continues.** These controls are not completely independent, because the added operators also bind BM3R1.

The distinction between new translation and accumulated protein also changed how we interpreted success. At baseline, late-LR stable protein remained at approximately **32.0% of peak**, despite new translation meeting the shutdown criterion. A future requirement for rapid protein removal would therefore need an output-specific degradation strategy.

The next engineering question was whether this behavior persisted beyond the nominal conditions, and which parameter relationships separated useful shutdown from failure.

## Engineering Cycle 2 — From Functional Shutdown to Defined Operating Constraints

### Design

**Design goal: identify the conditions under which delayed shutdown preserves a useful output window, and convert its failure modes into design constraints.**

We examined two sources of variation: growth-dependent changes in the coupled system and simultaneous variation of shutdown parameters. Because sRNA is consumed with its target, we also tested whether its effective transcriptional supply was sufficient relative to target-mRNA production.

### Build

For growth analysis, we regenerated oscillator input trajectories at **40-, 50-, and 60-minute doubling times**, updated protein dilution in the coupled model, and simulated each condition with its own no-shutdown control.

For parameter analysis, we implemented one-at-a-time sensitivity scans and a **512-point Latin hypercube sample across ten shutdown parameters**, combining literature-based ranges with explicitly declared design and context assumptions. Each sample was assessed against a no-shutdown reference matched to its target-mRNA transcription rate. This joint analysis used the nominal counter trajectory; full coupling and operator-load effects were assessed separately.

We additionally constructed a **225-point sRNA–target-mRNA transcription-supply scan**, using the same shutdown criteria as Cycle 1.

### Test

All three tested growth conditions met the shutdown and counter criteria. Translation-flux FWHM was **5.48, 6.98, and 8.54 h** at doubling times of **40, 50, and 60 min**, respectively. Thus, the tested growth changes preserved function but shifted absolute timing.

The broader joint parameter analysis gave a more conditional result: **277/512 samples (54.1%)** met the combined shutdown criteria, while **235 did not**. This demonstrates a feasible region within the declared parameter envelope, not universal robustness or a wet-lab success probability.

The transcription-supply scan helped explain this dependence. **124/225 combinations** passed. The lowest passing effective sRNA-to-target transcription ratio on this grid was **1.39**, and the median minimum passing ratio across target-mRNA levels was **2.68**. These values describe the tested grid, not universal biological thresholds.

In the joint parameter sample, **91.7%** of combinations with effective supply ratios between **4 and 8** passed. Raising the ratio further did not guarantee better performance: some combinations suppressed the initial peak too strongly. Sensitivity analysis also identified sRNA transcription strength as a major determinant of pulse width and late suppression.

![Self-consistent growth analysis](assets/05_self_consistent_growth.png)

*Each growth condition includes a regenerated oscillator trajectory and its own no-shutdown control.*

![Joint parameter analysis](assets/07_shutdown_joint_uncertainty.png)

*The Latin hypercube analysis measures coverage of a declared parameter envelope and reveals both passing and failing combinations.*

![Transcription-supply constraints](assets/04_transcription_balance_phase_map.png)

*The supply scan tests whether sRNA production can support sufficient target-mRNA removal while preserving the output peak.*

### Learn

The model supported retaining BM3R1–sRNA shutdown as a candidate architecture, with explicit operating constraints. Function was preserved across the three tested growth conditions, but only part of the broader parameter envelope supported both late shutdown and an adequate initial peak.

The main design consequence was to treat **effective sRNA-to-target transcriptional supply** as an explicit design variable. Ratios of **4–8** provide a model-supported starting range for future construct screening, while promoter strength and repression threshold offer ways to adjust timing. Increasing cassette dosage also increases BM3R1 operator load, so it cannot be treated as an isolated means of strengthening shutdown.

The next design cycle should therefore select and calibrate a specific sRNA–target construct, measure timing under the intended growth conditions, and define its connection to the target-count decoder. The present work establishes conditional feasibility for local shutdown of new translation; it does not yet demonstrate a complete count-to-N output system or experimental shutdown.

---

## Internal evidence record — not Wiki body text

This draft uses the current v3.3 delivery and does not introduce new simulations or experimentally validated claims.

- [Current model documentation](../README.md)
- [Delivery manifest and scope](../DELIVERY_MANIFEST.md)
- [Core results](../outputs/shutdown_analysis_20260909_150254/RESULTS.md)
- [Extended results](../outputs/shutdown_extended_20260909_133421/EXTENDED_RESULTS.md)
- [Growth-condition results](../outputs/shutdown_extended_20260909_133421/05_self_consistent_growth.csv)
- [Declared parameter ranges](../outputs/shutdown_extended_20260909_133421/07_parameter_ranges.csv)
- [Supply-ratio pass rates](../outputs/shutdown_extended_20260909_133421/07_supply_ratio_pass_rate.csv)
- [Existing Wiki model report, parameter provenance, and references](BM3R1驱动的sRNA_Shutdown模块_Wiki版报告.md)

Editorial decisions: omitted the proposed historical claim about discarded DNA-intermediate designs because the reviewed shutdown records do not establish that comparison; replaced blanket robustness claims with the actual pass fraction; distinguished translation suppression from protein clearance; retained the single-bit scope; and described the two cycles as an organization of the completed dry-lab work rather than asserting an undocumented chronology.

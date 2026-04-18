# qPCR Analysis Shiny App — Design Spec

**Date:** 2026-04-18  
**Status:** Approved

---

## Overview

A Shiny app for qPCR relative expression analysis, deployed to **GitHub Pages for free** using **Shinylive** (R runs entirely in the browser via WebAssembly — no server required). The user manually types or pastes Ct values into an editable table, selects analysis options, and gets a bar plot with statistical annotations plus downloadable outputs (PNG, PDF, CSV).

---

## Architecture

Single-file Shiny app (`app.R`) with:
- `ui` — sidebar layout (inputs left, outputs right)
- `server` — reactive ΔΔCt calculation, plot rendering, stat testing, file downloads

**Deployment:** The `shinylive` R package exports the app as a static site (HTML/JS/WASM). The exported `docs/` folder is committed to GitHub and served via GitHub Pages — free, no server, no Shinyapps.io account needed.

No external database. All state lives in the session (reactive values).

**R packages required:** `shiny`, `ggplot2`, `rstatix`, `DT`, `ggpubr`, `shinylive` (build-time only)

---

## Section 1 — Layout & Data Entry

**Sidebar inputs:**
- Reference gene name — text input (e.g. `GAPDH`)
- Control group name — text input (e.g. `Control`)
- Error bar type — radio buttons: `SD` or `SEM`
- "Paired samples" checkbox — unchecked by default
- Editable data table (via `DT` with editor) with columns:
  - `Sample` — sample identifier (used as pairing key when paired is checked)
  - `Group` — experimental group label
  - `Gene` — target gene name
  - `Ct_target` — Ct value for the gene of interest
  - `Ct_reference` — Ct value for the housekeeping/reference gene
- "Add Row" button — appends a blank row to the table
- "Clear Table" button — resets the table to empty

**Main panel:**
- Plot output area (one plot per gene, or faceted if multiple genes)
- Statistics results table
- Download buttons

---

## Section 2 — Analysis & Statistics

**ΔΔCt calculation (fully automatic, triggered reactively on table edits):**

1. ΔCt = `Ct_target` − `Ct_reference` (per sample)
2. ΔΔCt = ΔCt_sample − mean(ΔCt of control group samples)
3. Fold change = 2^(−ΔΔCt)
4. Control group samples always yield fold change = 1.0

**Statistical testing (auto-detected by group count, with pairing option):**

- When "Paired samples" is checked and there are 2 groups, a **paired t-test** is used. Unpaired is the default.
- For 3+ groups, pairing is not supported — ANOVA is always used regardless of the checkbox.
- The `Sample` column is used as the pairing identifier — each sample name must appear exactly once per group for pairing to work correctly.

| Groups | Paired? | Test | Post-hoc |
|--------|---------|------|----------|
| 2 | No | Unpaired t-test on ΔCt values | — |
| 2 | Yes | Paired t-test on ΔCt values | — |
| 3+ | — | One-way ANOVA on ΔCt values | Tukey HSD |

**Results table columns:**
- `Gene`, `Comparison` (e.g. `Treatment vs Control`), `p-value` (exact numeric, e.g. `0.023`), `Significant` (Yes / No, threshold p < 0.05)
- Significant rows highlighted in the table
- No significance stars — actual p-values only

**Multi-gene support:**
- If the data table contains more than one unique `Gene` value, analysis and plotting are run independently per gene.

---

## Section 3 — Plot

**Type:** Bar plot (one per gene, or faceted)

**Aesthetics:**
- X-axis: Group labels
- Y-axis: Mean fold change (2^−ΔΔCt)
- Error bars: mean ± SD or mean ± SEM (user-selected)
- Control group bar fixed at 1.0
- Bracket annotations above significant comparisons, labelled with exact p-value (e.g. `p = 0.023`)

**Color palette: Okabe-Ito (color-blind friendly)**

The Okabe-Ito palette is the standard for scientific publications and is safe for all common forms of color vision deficiency (deuteranopia, protanopia, tritanopia):

```r
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
```

Groups are assigned colors in order from this palette. The palette is applied via `scale_fill_manual()`.

**Rendered with:** `ggplot2` + `ggpubr` for bracket annotations

---

## Section 4 — Downloads

| Button | Output | Format |
|--------|--------|--------|
| Download Plot (PNG) | Current plot | PNG, 300 dpi |
| Download Plot (PDF) | Current plot | Vector PDF |
| Download Stats (CSV) | Statistics results table | CSV |

---

## Section 5 — Deployment (GitHub Pages via Shinylive)

1. Write `app.R` in the repo root.
2. Run `shinylive::export(".", "docs")` — this converts the app to static files in `docs/`.
3. Commit and push `docs/` to GitHub.
4. In the repo Settings → Pages, set source to `docs/` folder on `main` branch.
5. App is live at `https://<username>.github.io/<repo-name>/` — free, no server.

Re-export and push `docs/` whenever `app.R` changes.

---

## Error Handling

- If required columns are empty or non-numeric Ct values are entered, the app shows an inline error message and suppresses plot/stats output.
- If the control group name typed in the sidebar does not match any `Group` value in the table, an informative warning is shown.
- If only one sample per group is present, t-test/ANOVA cannot run — the app shows a message: "At least 2 samples per group required for statistical testing."
- If "Paired samples" is checked but sample names don't match across groups, a warning is shown and the test falls back to unpaired.

---

## Out of Scope

- File upload (CSV/Excel) — not needed; data entry is manual
- Multi-factor / two-way ANOVA — not in initial version
- User authentication or data persistence between sessions

Nice — this is *exactly* the right first move: **gather everything into one place + write an inventory**. That’s the hard psychological barrier, and you cleared it.

Now we want to turn what you wrote into a **migration plan with “thin wrappers”** so you can start getting wins without designing the whole universe.

Below is a roadmap doc you can paste into `ROADMAP.md` (or keep in `docs/roadmap.md`). It’s intentionally minimal and migration-focused.

---

## Gecko Roadmap (Migration-First, Minimal Steps)

### Mission for v0

Unify legacy MADNESS + Dalton parsing + your unit-sphere analysis tooling under a single namespace (`gecko`) **without rewriting**. The first goal is “boring and predictable loading.”

### Non-goals (for now)

* No agents that submit jobs.
* No generalized multi-step workflow engine.
* No new database system.
* No Trame rewrite (just relocate + make it runnable).

---

# Phase 0 — Inventory + boundaries (done + 1 small add)

✅ You created `migration/` and wrote `inventory.md`. Great.

**Add one more thing**: a simple “owner map” so you don’t mix concerns:

* **Parsing**: reading outputs → Python objects
* **Compute**: unit-sphere mapping + error metrics
* **Viz/apps**: plots + Trame UI
* **Workflow**: Slurm templates, step orchestration (defer)

Deliverable:

* `migration/inventory.md` includes a “belongs to” label per folder/file.

---

# Phase 1 — Create the gecko skeleton (tiny, 30–60 min)

Create:

```
src/gecko/
  __init__.py
  core/
    model.py
    load.py
  plugins/
    madness/
      __init__.py
      loader.py
    dalton/
      __init__.py
      loader.py
```

And a single public API:

* `gecko.load_calc(path)` (even if it returns a very simple object at first)

**Definition of success**

* `import gecko`
* `gecko.load_calc("/some/run")` returns something printable

---

# Phase 2 — Migrate parsing first (your #1 value)

You already said it clearly:

> “Ultimately, we only want a single dalton.py and madness.py for parsing.”

So let’s do it in a way that doesn’t break your working code.

## 2.1 Adopt a “thin wrapper” rule

* Keep your legacy parsers as-is in `migration/…`
* Copy them into `src/gecko/plugins/<code>/legacy/`
* Add one new wrapper module per code:

  * `src/gecko/plugins/madness/parse.py`
  * `src/gecko/plugins/dalton/parse.py`

These wrappers will:

* call legacy parsers
* normalize the result into a **minimal** common structure

## 2.2 Minimal common return type (don’t overdesign)

For migration, you only need:

```python
Calculation(
  code="madness" | "dalton",
  root=Path,
  artifacts={...},   # discovered files
  data={...},        # parsed content (raw/legacy allowed initially)
  meta={...}         # convenience info
)
```

**Definition of success**

* You can load 3 MADNESS runs and 3 Dalton runs using the same `gecko.load_calc`.

---

# Phase 3 — Standardize the “typical workflow” (but only as a recipe)

You pointed to `make_csv_data.ipynb` as the workflow.

We won’t port the notebook yet. Instead we make a **scriptable recipe**:

### Add `src/gecko/recipes/shg_csv.py`

A single function:

* `build_shg_dataframe(iter_runs(...)) -> pd.DataFrame`

Backed by:

* `iter_runs(root, code="madness|dalton", pattern=...)`
* `parse_beta(calc)` (calls your wrappers)

**Definition of success**

* One command or one Python snippet can rebuild your `shg_ijk.csv` from a directory tree.

This gives you a stable foundation for later GUIs/agents without committing to a database system.

---

# Phase 4 — Migrate `viz` (keep your best tool intact)

Your unit-sphere work is valuable *now*, so we keep it working.

## 4.1 Move without refactor

Copy:

* `migration.viz/field_error.py` → `src/gecko/viz/metrics.py`
* your unit-sphere mapping utilities → `src/gecko/viz/unit_sphere.py`
* trame `application.py` → `src/gecko/apps/trame_shg_viewer.py`

## 4.2 Replace direct CSV assumptions (later)

For now, let the Trame app still take `shg_ijk.csv`.

Later (not now), you’ll allow it to accept a `DataFrame` produced by your recipe.

**Definition of success**

* You can run the Trame app from the new repo with no behavioral changes.

---

# Phase 5 — Workflow tools: keep, but quarantine

Your `migration.cli` and `migration.db` are useful, but they are the biggest scope trap.

Do this:

* Copy them into `src/gecko/workflows/legacy/` (or keep them in migration for now)
* Add *only one* thing initially:

  * a tiny helper to **render** a Dalton input from template + molecule geometry
  * no submission, no polling, no job management (yet)

**Definition of success**

* You can reuse your templates cleanly from gecko code, but you’re not “building an HPC platform” yet.

---

# Phase 6 — Agent layer (placeholder only)

You absolutely should have an `agents/` folder, but in migration-v0 it contains only:

* `gecko/agents/README.md` describing the future tool interfaces
* `WorkflowSpec` dataclass (optional)

No execution. No LLM integration. Just a plan.

This prevents “agent gravity” from derailing migration.

---

## Immediate next steps (do these now)

### Step A — Create two “golden test fixtures”

Pick:

* 1 MADNESS run directory that you trust
* 1 Dalton run directory that you trust

Add:

* `tests/fixtures/madness_run/` (can be minimal subset)
* `tests/fixtures/dalton_run/`

Even if you can’t include full outputs, include the smallest set of files your parser needs.

### Step B — Implement `can_load()` detection

For each:

* MADNESS detection: e.g. presence of `responses/metadata.json` or your calcinfo JSON
* DALTON detection: presence of `.out` with signature lines, or your parsed json output

### Step C — Wrap legacy parsers

Your wrappers should do *only*:

* locate files
* call legacy parsing functions
* return the `Calculation` object

---

## Two questions that will save you time (but you can answer briefly)

1. For MADNESS, what is the single most reliable “marker file” in your runs?
   Examples: `responses/metadata.json`, `calcinfo.json`, `madness.json`, etc.

2. For Dalton, what’s the single output artifact you consider canonical?
   Examples: `DALTON.OUT`, `*.out`, `SIRIUS.RSP`, an existing JSON you generate, etc.

If you answer those, I’ll write you a super concrete “Phase 1–2 checklist” with exact filenames and where to put the wrapper functions, *based on your real layouts*.

Also: if you want, paste a `tree -L 2 migration/` output (just filenames) and I’ll propose exactly how to map each chunk into `gecko/plugins/*`, `gecko/viz`, `gecko/apps`, and what to leave in `workflows/legacy` for now.

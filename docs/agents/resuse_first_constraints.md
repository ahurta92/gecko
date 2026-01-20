# Reuse-First Constraints (Gecko)

## Mandatory Reuse Rules
- You MUST route all SHG-related work through existing public entrypoints:
  - `gecko.recipes.shg_csv` (public functions such as `build_beta_table(...)` or equivalent)
  - `gecko.core.iterators` (calculation discovery/iteration utilities)
- You MUST NOT reimplement logic already present in these modules.
- You MUST NOT create parallel pipelines (e.g., new SHG builders in other packages) when an existing path is viable.

## Allowed Edit Surface
- Prefer edits inside:
  - `gecko.recipes.shg_csv`
  - `gecko.core.iterators`
  - their immediate helper functions
- Creating any new module requires explicit user approval before proceeding.

## Process Enforcement
- Required “Inventory & Plan” step before coding:
  - List the exact functions to reuse and their locations.
  - Describe the minimal changes required.
- Required “Minimal diff” policy:
  - Max 3 files touched.
  - Max 120 lines added in total.

## Testing Enforcement
- Tests must exercise the existing public SHG entrypoint (e.g., `build_beta_table(...)`), not an alternative path.
- At least one test must prove the SHG path uses iterator + parser integration.

## Proof-of-Reuse (Required in final response)
- Include the import lines used.
- Provide a short call-trace of existing functions invoked.
- Explicitly confirm that no new SHG pipeline was created.

## Escalation Clause (Optional)
- If a required hook is missing, propose the smallest extension **inside existing modules** and request approval before any larger refactor or new module creation.

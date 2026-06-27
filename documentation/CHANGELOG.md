# Changelog

Chronological log of meaningful changes to this project. **Read the top entry at the start
of a session** to catch up on what your collaborator did since you last worked.

This doc's job — and how it differs from the others:

| Doc | Job |
|---|---|
| **`CHANGELOG.md`** (this file) | *When / who / what / why* — chronological, append-only. |
| `water_landing_cbf_design.md` | *Current state* — the living theory/approach/assumptions/decisions. |
| `water_landing_cbf_math.md` | *The math* — formal derivations. |
| `archive/` | Superseded historical docs. |

## Conventions

- **Newest entry on top.** Append above the existing entries; never rewrite past ones.
- One entry per working session (or per logically-grouped change). Keep it skimmable.
- Use the template below. Date is `YYYY-MM-DD`; author is a first name so collaborators can
  tell each other apart.

```
## YYYY-MM-DD — <author> — <short title>
**Branch/commit:** <branch> (<short-sha> if committed)
**What changed:** <the concrete change>
**Why:** <motivation / what problem it solves>
**Follow-ups / notes for collaborator:** <anything they should know or pick up>
**Files touched:** <key paths>
```

### Collaborator setup (one-time, per developer)

Each developer keeps their own local `CLAUDE.md` (gitignored, not shared via git). Add this
line to it so your Claude session reads this changelog at startup:

> At the start of each session, read `documentation/CHANGELOG.md` to see what changed since
> the last session before doing any work.

---

## 2026-06-27 — Jack — Set up `corbin-dev` branch + collaboration docs

**Branch/commit:** corbin-dev (off `lon-cbf-water-landing` @ `2b610a3`)
**What changed:**
- Created `corbin-dev` from the latest `lon-cbf-water-landing` (the augmented-longitudinal
  HOCBF water-landing module with exact autodiff/Taylor-jet Lie derivatives).
- Gitignored personal Claude context (`CLAUDE.md`, `PROJECT_CONTEXT.md`).
- Added this changelog and a living design doc (`water_landing_cbf_design.md`) that supersedes
  the dated implementation notes.
- Moved historical docs (`water_landing_cbf.md`, `water_landing_cbf_implementation_notes.md`)
  into `documentation/archive/`.
**Why:** Establish a shared working branch and a lightweight way for two parallel
Claude-using collaborators to stay in sync, plus a single maintained source of truth for the
design instead of a point-in-time session snapshot.
**Follow-ups / notes for collaborator:**
- Add the one-line changelog pointer to your own local `CLAUDE.md` (see Collaborator setup).
- Open design questions captured in `water_landing_cbf_design.md`: force-based `a_brk(V)`,
  the §3.3 contact-force barrier, and an independent `L_f³b` drift check.
**Files touched:** `.gitignore`, `documentation/CHANGELOG.md`,
`documentation/water_landing_cbf_design.md`, `documentation/archive/`.

# ADR-003: MIT licence

- **Status**: Accepted
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no

## Context

The repository contradicted itself: the README declared MIT in a badge and in its License section,
while the `LICENSE` file contained the GNU General Public License version 2, which is what GitHub
reported. A user could not tell which terms applied.

## Decision

rpsHunter is released under the MIT licence. The `LICENSE` file is replaced with the MIT text and a
`CITATION.cff` is added. Every commit in the history up to this change was authored by the copyright
holder, so no other contributor's agreement was needed.

## Consequences

- Positive: one unambiguous, permissive licence, the same as the sibling pipeline RetroSeek.
- Negative: derived works are no longer required to remain open source.
- Neutral: copies obtained while the GPL v2 text was in place were offered under those terms and
  remain so; MIT applies from this change onward.

## Alternatives considered

- **Keep GPL v2 and correct the README**: preserves copyleft; not chosen because the stated
  intent, in the README from the start, was MIT.
- **GPL v3 or later**: copyleft in its current form; same reason.

## Revisit trigger

A contribution arrives under terms that are incompatible with MIT, or a dependency's licence
requires a different one.

## References

- `LICENSE`, `CITATION.cff`, `README.md`.

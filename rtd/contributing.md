# Contributing

This page gives initial guidance for contributing to ICARION.

## Development philosophy

ICARION prioritizes:

- physical correctness,
- reproducibility,
- modularity,
- clear configuration,
- [validation](validation.md) before feature expansion,
- and transparent output formats.

## Suggested workflow

1. Open or discuss an issue describing the proposed change.
2. Create a feature branch.
3. Add or update tests where appropriate.
4. Run the relevant validation or regression checks.
5. Open a pull request with a clear summary of the change.

## Documentation changes

Documentation changes should be made in the `rtd/` directory for Read the Docs pages. Technical Markdown files elsewhere in the repository may be used for lower-level reference notes, but the public user manual should remain readable and beginner-friendly.

## Adding new physical models

When adding a new physical model, include:

- a clear description of the assumptions,
- configuration fields and defaults,
- update allows configuration schema, if applicable,
- unit conventions,
- at least one small test case,
- and one validation or sanity-check example when possible.

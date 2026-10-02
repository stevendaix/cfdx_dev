# CFDX Documentation

CFDX documentation is organized as four complementary domains:

| Domain | Question | Primary audience |
|---|---|---|
| [Theory](theory/README.md) | Why does the method work, and what is the mathematics? | CFD/numerical-method users and developers |
| [User](user/README.md) | How do I use CFDX? | CFD users |
| [Developer](developer/README.md) | How is CFDX implemented and extended? | Developers |
| [V&V](vv/README.md) | What evidence demonstrates correctness and qualification? | Developers, reviewers, users |

The domains are deliberately separated. A page must have one primary role and one authoritative owner.

## Documentation flow

```text
THEORY
  mathematical definition and numerical properties
        |
        v
DEVELOPER
  implementation contract and software architecture
        |
        v
V&V
  executable verification, validation and qualification evidence
        |
        v
USER
  supported capability and practical usage
```

A link may point upstream or downstream, but the same technical fact must not be independently maintained in multiple domains.

## Repository structure

```text
docs/
├── README.md
├── theory/                 # executable scientific course
├── user/                   # task-oriented user documentation
├── developer/              # architecture and implementation documentation
├── vv/                     # V&V portal and methodology
├── validation/             # current V&V evidence store; migrated incrementally
└── references/
    └── bibliography.bib    # authoritative bibliography
```

The existing `docs/validation/` tree remains in place during migration. It is not duplicated or mass-renamed in the foundation phase.

## Documentation principles

1. **One fact, one owner.**
2. **Theory is explanatory and mathematical, not a user manual.**
3. **User documentation is task-oriented, not an implementation guide.**
4. **Developer documentation describes contracts and architecture, not numerical theory already owned by Theory.**
5. **V&V documentation is evidence-driven and never promotes status without retained evidence.**
6. **Executable scientific content is preferred over static numerical claims.**
7. **Generated figures and tables have a reproducible source.**
8. **Bibliographic metadata has one authoritative source.**
9. **Historical documents are clearly identified as historical and cannot silently become current status authorities.**
10. **Documentation structure is itself tested by CI.**

## Planned publication model

The source documentation is intended to support:

- an interactive web documentation site;
- executable Python scientific examples;
- interactive visualisations where useful;
- generated V&V figures and tables;
- Typst-generated PDF books and technical reports.

Markdown remains useful for prose, but it is not the scientific execution layer. Python is the preferred executable layer for theory experiments and reproducible numerical demonstrations.

## Migration policy

The migration is incremental:

1. establish ownership and navigation;
2. preserve existing evidence;
3. introduce the publication/tooling layer;
4. migrate one documentation family at a time;
5. validate links and generated content;
6. remove obsolete duplicates only after their replacement is proven.

Do not perform a destructive documentation rewrite in a single PR.

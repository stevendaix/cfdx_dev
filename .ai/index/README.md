# CFDX Code Intelligence Index

The planned code-intelligence layer provides structural repository queries to AI agents.

Target queries include symbol definitions, references, callers/callees, include/dependency relationships, inheritance, module impact, tests covering implementation and validation cases exercising implementation.

Candidate technologies include Clang AST/libclang, Python AST, Tree-sitter, compile_commands.json, LSP information and ripgrep. A local SQLite index is a candidate implementation.

The index must be reproducible and incrementally updateable and must never become the authoritative source of repository truth.

Text search answers where text occurs. Code intelligence should additionally answer structural questions such as what calls a function or which implementation depends on a class.
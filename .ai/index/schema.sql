-- CFDX code-intelligence logical schema.
-- Generated indexes are disposable derived data; this schema is not runtime CFDX state.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS index_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS files (
    id INTEGER PRIMARY KEY,
    path TEXT NOT NULL UNIQUE,
    language TEXT NOT NULL,
    content_hash TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    revision TEXT,
    indexed_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS symbols (
    id INTEGER PRIMARY KEY,
    file_id INTEGER NOT NULL REFERENCES files(id) ON DELETE CASCADE,
    kind TEXT NOT NULL,
    qualified_name TEXT NOT NULL,
    name TEXT NOT NULL,
    line_start INTEGER NOT NULL,
    column_start INTEGER,
    line_end INTEGER,
    column_end INTEGER
);

CREATE INDEX IF NOT EXISTS idx_symbols_file ON symbols(file_id);
CREATE INDEX IF NOT EXISTS idx_symbols_qualified_name ON symbols(qualified_name);

CREATE TABLE IF NOT EXISTS code_references (
    id INTEGER PRIMARY KEY,
    source_symbol_id INTEGER REFERENCES symbols(id) ON DELETE CASCADE,
    target_symbol_id INTEGER REFERENCES symbols(id) ON DELETE CASCADE,
    file_id INTEGER NOT NULL REFERENCES files(id) ON DELETE CASCADE,
    kind TEXT NOT NULL,
    line INTEGER NOT NULL,
    column INTEGER
);

CREATE INDEX IF NOT EXISTS idx_code_references_source ON code_references(source_symbol_id);
CREATE INDEX IF NOT EXISTS idx_code_references_target ON code_references(target_symbol_id);
CREATE INDEX IF NOT EXISTS idx_code_references_file ON code_references(file_id);

CREATE TABLE IF NOT EXISTS dependencies (
    id INTEGER PRIMARY KEY,
    source_file_id INTEGER NOT NULL REFERENCES files(id) ON DELETE CASCADE,
    target_path TEXT NOT NULL,
    kind TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_dependencies_source ON dependencies(source_file_id);
CREATE INDEX IF NOT EXISTS idx_dependencies_target ON dependencies(target_path);

CREATE TABLE IF NOT EXISTS tests (
    id INTEGER PRIMARY KEY,
    file_id INTEGER NOT NULL REFERENCES files(id) ON DELETE CASCADE,
    name TEXT NOT NULL,
    framework TEXT,
    line_start INTEGER,
    line_end INTEGER
);

CREATE INDEX IF NOT EXISTS idx_tests_file ON tests(file_id);

CREATE TABLE IF NOT EXISTS validation_cases (
    id INTEGER PRIMARY KEY,
    path TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL,
    category TEXT,
    status TEXT,
    source_revision TEXT
);

CREATE TABLE IF NOT EXISTS test_links (
    test_id INTEGER NOT NULL REFERENCES tests(id) ON DELETE CASCADE,
    symbol_id INTEGER NOT NULL REFERENCES symbols(id) ON DELETE CASCADE,
    relation TEXT NOT NULL,
    PRIMARY KEY (test_id, symbol_id, relation)
);

CREATE TABLE IF NOT EXISTS validation_links (
    validation_case_id INTEGER NOT NULL REFERENCES validation_cases(id) ON DELETE CASCADE,
    symbol_id INTEGER NOT NULL REFERENCES symbols(id) ON DELETE CASCADE,
    relation TEXT NOT NULL,
    PRIMARY KEY (validation_case_id, symbol_id, relation)
);

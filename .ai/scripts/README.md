# CFDX AI Scripts

Planned maintenance utilities:

- validate_skill.py — validate skill metadata and structure.
- validate_ai_config.py — validate AI framework configuration.
- index_repository.py — initial repository index generation.
- update_index.py — incremental index update.
- generate_agent_configs.py — generate vendor adapters from canonical definitions.

Scripts must fail explicitly on malformed configuration and must not silently rewrite project files.
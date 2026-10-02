# %% [markdown]
# # Qualification matrix
#
# This executable schema contains no fabricated results. It converts retained campaign evidence into bounded qualification records.
#
# %%
from dataclasses import dataclass
from typing import Literal
#
# %%
Status = Literal["IMPLEMENTED", "VERIFIED", "VALIDATED", "QUALIFIED",
                  "FAILED", "NOT_QUALIFIED", "NOT_RUN"]
#
@dataclass(frozen=True)
class EvidenceRow:
    capability: str
    population: str
    requirement: str
    campaign: str
    metric: str
    criterion: str
    evidence_artifact: str
    status: Status
#
def qualifies(row: EvidenceRow) -> bool:
    return bool(row.evidence_artifact) and row.status == "QUALIFIED"

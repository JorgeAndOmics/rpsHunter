"""The pipeline's code knows no domain family by name.

Families are defined in the config by accession (ADR-012). A family name written
into a script ties the tool to one protein family and breaks when a database
renames a model, so none may appear in the workflow's code, comments included.
The names checked are those of the shipped example and of its models.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FAMILY_NAME = re.compile(r"KRAB|SSXRD|zf-C2H2|zf-H2C2|ZnF_C2H2|ZMS1|PRDM|\bSET\b")


def test_no_family_name_appears_in_the_workflow_code():
    scripts = REPO / "workflow" / "scripts"
    sources = [REPO / "rpsHunter", REPO / "workflow" / "Snakefile"]
    for pattern in ("*.py", "*.R", "*.sh"):
        sources += sorted(scripts.rglob(pattern))
    offenders = [
        f"{path.relative_to(REPO)}:{number}: {line.strip()[:80]}"
        for path in sources
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if FAMILY_NAME.search(line)
    ]

    assert offenders == []

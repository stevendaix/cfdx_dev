from cfdx import CFDXSession
from cfdx.application import Application
from cfdx.application.commands import ConfigureProbes
from cfdx.probe import Probe, ProbeCatalog, ProbeField


def test_configure_probes_command_attaches_probes_to_case() -> None:
    session = CFDXSession()
    application = Application(session)

    assert session.case.probes == []
    catalog = ProbeCatalog([Probe("p_wall", 1.0, 2.0, 3.0, ProbeField.PRESSURE)])

    state = application.execute(ConfigureProbes(catalog))

    assert session.case.probes == list(catalog.probes)
    assert session.case_revision == 1
    assert state.case_revision == 1
    assert not session.requires_restart


def test_configure_probes_replaces_existing_catalogue() -> None:
    session = CFDXSession()
    application = Application(session)

    application.execute(ConfigureProbes(ProbeCatalog([
        Probe("a", 0, 0, 0, ProbeField.PRESSURE),
        Probe("b", 1, 1, 1, ProbeField.U_X),
    ])))
    assert len(session.case.probes) == 2

    application.execute(ConfigureProbes(ProbeCatalog([Probe("only", 2, 2, 2, ProbeField.U_MAGNITUDE)])))
    assert [p.name for p in session.case.probes] == ["only"]


def test_configure_probes_does_not_require_solver_restart() -> None:
    session = CFDXSession()
    application = Application(session)
    application.execute(ConfigureProbes(ProbeCatalog([Probe("p", 0, 0, 0)])))
    assert not session.requires_restart


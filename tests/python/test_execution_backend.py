from cfdx.execution_backend import SlurmConfig, slurm_command


def test_slurm_command_is_explicit() -> None:
    command = slurm_command(["cfdx_solver", "case.json"], SlurmConfig(partition="cfd", tasks=8))
    assert command[:3] == ("sbatch", "--wait", "--nodes")
    assert "--partition" in command
    assert "cfdx_solver case.json" in command


def test_empty_slurm_command_rejected() -> None:
    try:
        slurm_command([], SlurmConfig())
    except ValueError:
        pass
    else:
        raise AssertionError("empty command must be rejected")

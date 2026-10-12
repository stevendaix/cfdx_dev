from cfdx.output import OutputThrottle


def test_throttle_batches_and_preserves_stream() -> None:
    emitted = []
    throttle = OutputThrottle(lambda line, is_stderr: emitted.append((line, is_stderr)), interval=0.05, max_lines=2)
    throttle.push("a", False)
    throttle.push("b", True)
    assert throttle.flush(now=1.0) == 2
    assert emitted == [("a", False), ("b", True)]
    throttle.push("c", False)
    assert throttle.flush(now=1.01) == 0
    assert throttle.pending() == 1
    assert throttle.flush(now=1.06) == 1
    assert emitted[-1] == ("c", False)


def test_throttle_validates_limits() -> None:
    emitted = []
    for interval, max_lines in ((0, 10), (0.1, 0)):
        try:
            OutputThrottle(emitted.append, interval=interval, max_lines=max_lines)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid limits must raise ValueError")

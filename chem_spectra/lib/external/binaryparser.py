from typing import Callable, Optional, Tuple


def get_openlab_readers() -> Tuple[Optional[Callable], Optional[Callable]]:
    try:
        try:
            from BinaryParser import read_lc  # type: ignore
        except Exception:
            from BinaryParser.openlab import read_lc  # type: ignore
        try:
            from BinaryParser.openlab import read_ms  # type: ignore
        except Exception:
            from BinaryParser import read_ms  # type: ignore
        return read_lc, read_ms
    except Exception:
        return None, None


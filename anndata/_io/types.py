from typing import AnyStr

try:
    from typing import Protocol
except ImportError:
    from typing_extensions import Protocol


class SupportsRead(Protocol[AnyStr]):
    def read(self, len: int) -> AnyStr:
        ...


class SupportsWrite(Protocol[AnyStr]):
    def write(self, content: AnyStr) -> int:
        ...

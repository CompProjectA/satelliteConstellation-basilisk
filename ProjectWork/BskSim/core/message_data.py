# message_data.py
from dataclasses import dataclass
from typing import Any, Union

@dataclass(frozen=True)
class MessageData:
    """
    Lightweight container for inter-satellite messages.

    Fields:
        message_content: The message payload (text).
        timeSent:        Simulation time when the message was sent (minutes, or any numeric).
        objectActive:    The target/peer object involved in the message (e.g., leader/child instance).

    Equality:
        Two MessageData objects are equal only if all three fields match,
        including the exact same 'objectActive' reference.
    """
    message_content: str
    timeSent: Union[float, int]
    objectActive: Any

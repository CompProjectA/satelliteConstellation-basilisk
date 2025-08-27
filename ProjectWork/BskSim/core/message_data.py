# message_data.py
class MessageData:
    def __init__(self, message_content: str, timeSent, objectActive):
        self.message_content = message_content
        self.timeSent = timeSent
        self.objectActive = objectActive

    def __eq__(self, other):
        if isinstance(other, MessageData):
            return (
                self.message_content == other.message_content
                and self.timeSent == other.timeSent
                and self.objectActive == other.objectActive
            )
        return False

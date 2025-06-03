
class MessageData:
    def __init__(self, message_content:str, timeSent, objectActive):
        
        self.message_content = message_content
        self.timeSent = timeSent
        self.objectActive = objectActive
        #Sending is sent by this object, Receiving means received from another object
    def __eq__(self, value):
        if isinstance(value,MessageData):
            return(
                self.message_content == value.message_content and
                self.timeSent == value.timeSent and
                self.objectActive == value.objectActive
            )
        return False
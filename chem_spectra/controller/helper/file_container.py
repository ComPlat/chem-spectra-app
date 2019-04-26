class FileContainer:
    def __init__(self, src=False):
        self.name = src and src.filename
        self.bcore = src and src.stream.read()
        self.core = self.bcore and self.bcore.decode('utf-8', errors='ignore')


    def from_str(self, core):
        self.core = core
        self.bcore = str.encode(core)
        return self

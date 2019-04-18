class FileContainer:
    def __init__(self, src=False):
        self.name = src and src.filename
        self.core = src and src.stream.read().decode('utf-8', errors='ignore')

    def from_str(self, core):
        self.core = core
        return self

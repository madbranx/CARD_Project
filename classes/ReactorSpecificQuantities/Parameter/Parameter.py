import casadi as casADi

class Parameter:
    def __init__(self, log, name, value=None):
        self.log = log
        self.name  = name
        self.value = None
        if value is not None:
            self.setValue(value)

    def setValue(self, value):
        self.value = value

    def getName(self):
        return self.name

    def getValue(self):
        return self.value
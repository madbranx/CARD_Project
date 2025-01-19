import casadi as casADi

class Parameter:
    def __init__(self, log, name, value=None):
        self.log = log
        self.name  = name
        if value is not None:
            self.value = casADi.MX.sym(name, value)

    def getName(self):
        return self.name

    def getValue(self):
        return self.value

    def setValue(self, value):
        self.value = casADi.MX.sym(self.name, value)
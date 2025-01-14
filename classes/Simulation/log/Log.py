import datetime

class Log:
    def __init__(self, name):
        self.text = ""
        self.currentText = ""  # text of the current timestep
        self.name = name
        self.currentWarnings = []
        self.warnings = ""
        self.currentErrors = []
        self.errors = ""
        ###########
        self.addEntry("Starting log of simulation \"" + self.name + "\"", 0)
        self.updateLog()

    def addEntry(self, text, level):
        indent = level * "\t"
        text = "{}{}\n".format(indent, text)
        # print(text)
        self.currentText += text

    def addWarning(self, text, level=1):
        indent = level*"\t"
        warning = "{}WARNING: {}\n".format(indent, text)
        print(warning)
        self.currentText += warning
        self.currentWarnings.append(text)

    def addError(self, text, level=1):
        indent = level*"\t"
        error = "{}ERROR: {}\n".format(indent, text)
        # print(error)
        self.currentText += error
        self.currentErrors.append(text)

    def updateLog(self, time=None):
        if self.currentText != "":
            if time is None:
                text = self.currentText
                entry = "{}\n".format(text)
            else:
                timeString = str(datetime.timedelta(seconds=time))
                text = "events:\n" + self.currentText
                entry = "\t[{}] {}\n".format(timeString, text)
            # print(entry)
            self.currentText = ""
            self.text += entry
            self.__updateWarnings(time)
            self.__updateErrors(time)
            # if self.errors:
            #     self.export()

    def __updateWarnings(self, time=None):
        if self.currentWarnings:
            if time is None:
                self.warnings += "\n"
            else:
                timeString = str(datetime.timedelta(seconds=time))
                self.warnings += "\n" + timeString + "\n"
            for warning in self.currentWarnings:
                self.warnings += "\t" + warning + "\n"
            self.currentWarnings = []

    def __updateErrors(self, time=None):
        if self.currentErrors:
            if time is None:
                self.errors += "\n"
            else:
                timeString = str(datetime.timedelta(seconds=time))
                self.errors += "\n" + timeString + "\n"
            for error in self.currentErrors:
                self.errors += "\t" + error + "\n"
            self.currentErrors = []

    def export(self):
        file_path = "logs/log_" + self.name + ".txt"
        file = open(file_path, 'w')
        if self.warnings:
            file.write("WARNINGS:" + self.warnings + "\n")
        else:
            file.write("no warnings found\n")
        if self.errors:
            file.write("ERRORS:" + self.errors + "\n")
        else:
            file.write("no errors found\n")
        file.write("\n" + "#"*60 + "\n")
        file.write(self.text)
        print("log exported to " + file_path)
        if self.errors:
            raise RuntimeError("Error occurred in simulation. see log for details")

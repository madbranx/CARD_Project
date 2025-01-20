import numpy as np

class Discretization:
    EQUIDISTANT = 1
    ARRAY = 2
    FUNCTION = 3
    def __init__(self, log, num_volumes, kind=EQUIDISTANT, **kwargs):
        self.log = log
        self.num_faces = num_volumes + 1
        self.kind = kind
        self.kwargs = kwargs
        self.start = kwargs.get('start', 0.0)
        self.end = kwargs.get('end', 0.0)
        self.faces = None
        self.centroids = None
        self.differences = None
        self.__create()

    def __create(self):
        if self.kind == self.EQUIDISTANT:
            self.faces = self.__equidistant()
        elif self.kind == self.ARRAY:
            self.faces = self.__array()
        elif self.kind == self.FUNCTION:
            self.faces = self.__function()
        else:
            pass # add log entry

        self.centroids =  (self.faces[:-1] + self.faces[1:]) / 2
        self.differences = np.diff(self.centroids)

    def __equidistant(self):
        return np.linspace(self.start, self.end, self.num_faces)

    def __array(self):
        pass  # TODO

    def __function(self):
        pass  # TODO

    def get_faces(self):
        return self.faces

    def get_centroids(self):
        return self.centroids

    def get_difference(self):
        return self.differences

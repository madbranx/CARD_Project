import numpy as np

class Discretization:
    EQUIDISTANT = 1
    ARRAY = 2

    def __init__(self, log, num_volumes, kind=EQUIDISTANT, **kwargs):
        self.log = log
        self.num_faces = num_volumes + 1
        self.kind = kind
        self.kwargs = kwargs
        self.start = kwargs.get('start', 0.0)
        self.end = kwargs.get('end', 1.0)
        self.faces = None
        self.centroids = None
        self.differences = None
        self.__create()

    def __create(self):
        self.log.addEntry("type = " + str(self.kind), 3)
        self.log.addEntry("start = " + str(self.start), 3)
        self.log.addEntry("end = " + str(self.end), 3)
        self.log.addEntry("volumes = " + str(self.num_faces - 1), 3)

        if self.kind == self.EQUIDISTANT:
            self.faces = self.__equidistant()
        elif self.kind == self.ARRAY:
            self.faces = self.__array()
        else:
            self.log.addWarning("Invalid discretization kind. Equidistant Discretization will be used!", 2)
            self.faces = self.__equidistant()

        self.centroids = (self.faces[:-1] + self.faces[1:]) / 2
        self.differences = np.diff(self.faces)

    def __equidistant(self):
        return np.linspace(self.start, self.end, self.num_faces)

    def __array(self):
        ranges = self.kwargs.get('ranges', [])  # Expect [(start, end, step), ...]
        if not ranges:
            self.log.addWarning("No ranges provided for ARRAY discretization. Equidistant Discretization will be used!", 2)
            return self.__equidistant()

        total_factor = sum(factor for _, _, factor in ranges)

        faces = []
        remaining_volumes = self.num_faces - 1

        for i, (start, end, factor) in enumerate(ranges):
            # Allocate volumes proportionally based on the factor
            n_volumes_range = max(1, int(round((factor / total_factor) * remaining_volumes)))
            remaining_volumes -= n_volumes_range
            if i == len(ranges) - 1:  # Ensure all volumes are used by the last range
                n_volumes_range += remaining_volumes

            segment = np.linspace(start, end, n_volumes_range + 1)
            faces.extend(segment[:-1])

        faces.append(ranges[-1][1])  # Ensure the last end is included
        return np.array(faces)


    def get_faces(self):
        return self.faces

    def get_centroids(self):
        return self.centroids

    def get_differences(self):
        return self.differences

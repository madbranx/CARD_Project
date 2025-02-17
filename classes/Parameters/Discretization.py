import numpy as np

class Discretization:
    EQUIDISTANT = 1
    ARRAY = 2
    RELATIVE_ARRAY = 3

    def __init__(self, num_volumes, kind=EQUIDISTANT, **kwargs):
        self.num_volumes = num_volumes
        self.kind = kind
        self.kwargs = kwargs
        self.start = kwargs.get('start', 0.0)
        self.end = kwargs.get('end', 1.0)
        self.faces = None
        self.centroids = None
        self.differences_faces = None

        self.__create()

    def __create(self):
        if self.kind == self.EQUIDISTANT:
            self.faces = self.__equidistant()
        elif self.kind == self.ARRAY:
            self.faces = self.__array()
        elif self.kind == self.RELATIVE_ARRAY:
            self.faces = self.__relative_array()
        else:
            self.faces = self.__equidistant()

        self.centroids = (self.faces[:-1] + self.faces[1:]) / 2
        self.differences_faces = np.diff(self.faces)
        self.differences_centroids = np.diff(self.centroids)

    def __equidistant(self):
        return np.linspace(self.start, self.end, self.num_volumes+1)

    def __array(self):
        ranges = self.kwargs.get('ranges', [])  # Expect [(start, end, step), ...]

        total_factor = sum(factor for _, _, factor in ranges)

        faces = []
        remaining_volumes = self.num_volumes

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

    def __relative_array(self):
        ranges = self.kwargs.get('ranges', [])  # Expect [(start, end, step), ...]

        total_length = self.end - self.start
        total_weight = sum((factor * (end - start) / total_length) for start, end, factor in ranges)

        faces = []
        remaining_volumes = self.num_volumes

        for i, (start, end, factor) in enumerate(ranges):
            segment_length = end - start
            segment_weight = (factor * segment_length) / total_length
            n_volumes_range = max(1, int(round((segment_weight / total_weight) * self.num_volumes)))
            remaining_volumes -= n_volumes_range

            if i == len(ranges) - 1:  # Ensure last segment gets any remaining volumes
                n_volumes_range += remaining_volumes

            segment = np.linspace(start, end, n_volumes_range + 1)
            faces.extend(segment[:-1])

        faces.append(ranges[-1][1])
        return np.array(faces)


    def get_faces(self):
        return self.faces

    def get_centroids(self):
        return self.centroids

    def get_differences_faces(self):
        return self.differences_faces

    def get_differences_centroids(self):
        return self.differences_centroids
